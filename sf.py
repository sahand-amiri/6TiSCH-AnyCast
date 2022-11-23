from __future__ import absolute_import
# =========================== imports =========================================

from builtins import range
from builtins import object
import random
import sys
from abc import abstractmethod

import netaddr

import SimEngine
from . import MoteDefines as d
from . import sixp

# =========================== defines =========================================

# =========================== helpers =========================================

# =========================== body ============================================

class SchedulingFunction(object):
    def __new__(cls, mote):
        settings    = SimEngine.SimSettings.SimSettings()
        class_name  = u'SchedulingFunction{0}'.format(settings.sf_class)
        return getattr(sys.modules[__name__], class_name)(mote)

class SchedulingFunctionBase(object):

    SLOTFRAME_HANDLE = 0

    def __init__(self, mote):

        # store params
        self.mote            = mote

        # singletons (quicker access, instead of recreating every time)
        self.settings        = SimEngine.SimSettings.SimSettings()
        self.engine          = SimEngine.SimEngine.SimEngine()
        self.log             = SimEngine.SimLog.SimLog().log

    # ======================= public ==========================================

    # === admin

    @abstractmethod
    def start(self):
        """
        tells SF when should start working
        """
        raise NotImplementedError() # abstractmethod

    @abstractmethod
    def stop(self):
        '''
        tells SF when should stop working
        '''
        raise NotImplementedError() # abstractmethod

    # === indications from other layers

    @abstractmethod
    def indication_neighbor_added(self, neighbor_mac_addr):
        pass

    @abstractmethod
    def indication_tx_cell_elapsed(self, cell, sent_packet):
        """[from TSCH] just passed a dedicated TX cell. used=False means we didn't use it.

        """
        raise NotImplementedError() # abstractmethod

    @abstractmethod
    def indication_rx_cell_elapsed(self, cell, received_packet):
        raise NotImplementedError() # abstractmethod

    @abstractmethod
    def indication_parent_change(self, old_parent, new_parent):
        """
        [from RPL] decided to change parents.
        """
        raise NotImplementedError() # abstractmethod

    @abstractmethod
    def detect_schedule_inconsistency(self, peerMac):
        raise NotImplementedError() # abstractmethod

    @abstractmethod
    def recv_request(self, packet):
        raise NotImplementedError() # abstractmethod

    @abstractmethod
    def clear_to_send_EBs_DATA(self):
        raise NotImplementedError() # abstractmethod


class SchedulingFunctionSFNone(SchedulingFunctionBase):

    def __init__(self, mote):
        super(SchedulingFunctionSFNone, self).__init__(mote)

    def start(self):
        pass # do nothing

    def stop(self):
        pass # do nothing

    def indication_neighbor_added(self, neighbor_mac_addr):
        pass # do nothing

    def indication_tx_cell_elapsed(self, cell, sent_packet):
        pass # do nothing

    def indication_rx_cell_elapsed(self, cell, received_packet):
        pass # do nothing

    def indication_parent_change(self, old_parent, new_parent):
        pass # do nothing

    def detect_schedule_inconsistency(self, peerMac):
        pass # do nothing

    def recv_request(self, packet):
        pass # do nothing

    def clear_to_send_EBs_DATA(self):
        # always return True
        return True


class SchedulingFunctionMSF(SchedulingFunctionBase):

    SLOTFRAME_HANDLE_AUTONOMOUS_CELLS = 1
    SLOTFRAME_HANDLE_NEGOTIATED_CELLS = 2
    DEFAULT_CELL_LIST_LEN = 5
    MAX_RETRY = 3
    TX_CELL_OPT   = [d.CELLOPTION_TX]
    RX_CELL_OPT   = [d.CELLOPTION_RX]
    NUM_INITIAL_NEGOTIATED_TX_CELLS = 1
    NUM_INITIAL_NEGOTIATED_RX_CELLS = 0

    def __init__(self, mote):
        # initialize parent class
        super(SchedulingFunctionMSF, self).__init__(mote)

        # (additional) local variables
        self.num_tx_cells_elapsed = 0
        self.num_tx_cells_used    = 0
        self.tx_cell_utilization  = 0
        self.num_rx_cells_elapsed = 0
        self.num_rx_cells_used    = 0
        self.rx_cell_utilization  = 0
        self.locked_slots         = set([]) # slots in on-going ADD transactions
        self.retry_count          = {}      # indexed by MAC address
        
        if self.mote.id == 0:
            self.schedulingMatrix = {}
            self.packetDroppedPerLink = {}

    # ======================= public ==========================================

    # === admin

    def start(self):
        # install slotframes for MSF, which have the same length as
        # Slotframe 0
        slotframe_0 = self.mote.tsch.get_slotframe(0)
        self.mote.tsch.add_slotframe(
            slotframe_handle = self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS,
            length           = slotframe_0.length
        )
        self.mote.tsch.add_slotframe(
            slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS,
            length           = slotframe_0.length
        )

        # install our autonomous RX cell
        self.allocate_autonomous_rx_cell()

        if self.mote.dagRoot:
            # do nothing
            pass
        else:
            self._housekeeping_collision()

    def stop(self):
        # uninstall the slotframes entirely instead of removing the
        # cells there one by one
        self.mote.tsch.delete_slotframe(self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS)
        self.mote.tsch.delete_slotframe(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)

        if self.mote.dagRoot:
            # do nothing
            pass
        else:
            self.engine.removeFutureEvent(
                (self.mote.id, u'_housekeeping_collision')
            )

    # === indications from other layers

    def indication_neighbor_added(self, neighbor_mac_addr):
        pass

    def indication_tx_cell_elapsed(self, cell, sent_packet):
        preferred_parent = self.mote.rpl.getPreferredParent()
        if (
                preferred_parent
                and
                (cell.mac_addr == preferred_parent)
                and
                (cell.options == [d.CELLOPTION_TX])
            ):
            self._update_cell_counters(self.TX_CELL_OPT, bool(sent_packet))
            # adapt number of cells if necessary
            if d.MSF_MAX_NUMCELLS <= self.num_tx_cells_elapsed:
                tx_cell_utilization = (
                    self.num_tx_cells_used /
                    float(self.num_tx_cells_elapsed)
                )
                if tx_cell_utilization != self.tx_cell_utilization:
                    self.log(
                        SimEngine.SimLog.LOG_MSF_TX_CELL_UTILIZATION,
                        {
                            u'_mote_id'    : self.mote.id,
                            u'neighbor'    : preferred_parent,
                            u'value'       : u'{0}% -> {1}%'.format(
                                int(self.tx_cell_utilization * 100),
                                int(tx_cell_utilization * 100)
                            )
                        }
                    )
                    self.tx_cell_utilization = tx_cell_utilization
                self._adapt_to_traffic(preferred_parent, self.TX_CELL_OPT)
                self._reset_cell_counters(self.TX_CELL_OPT)

    def indication_rx_cell_elapsed(self, cell, received_packet):
        preferred_parent = self.mote.rpl.getPreferredParent()
        if not preferred_parent:
            # nothing to do
            return

        if (
                (cell.mac_addr == preferred_parent)
                and
                (cell.options == [d.CELLOPTION_RX])
            ):
            self._handle_rx_cell_elapsed_event(bool(received_packet))
        elif (
                (cell.mac_addr == None)
                and
                (
                    cell.slotframe.slotframe_handle ==
                    self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS
                )
            ):
            if (
                received_packet is None
                or
                received_packet[u'mac'][u'srcMac'] == preferred_parent
                ):
                if not self.get_negotiated_rx_cells(preferred_parent):
                    self._handle_rx_cell_elapsed_event(bool(received_packet))
                else:
                    # ignore this notification
                    pass
            elif (
                self.get_negotiated_rx_cells(received_packet[u'mac'][u'srcMac'])
                ):
                self._handle_rx_cell_elapsed_event(False)
                assert cell.options == [d.CELLOPTION_RX]
                # we received a packet on our autonomous RX cell, with the
                # source mote of which we have negotiated RX cells. The
                # source mote must have lost the negotaited RX cells, TX
                # on its viewpoint. Remove them now.
                self._clear_cells(received_packet[u'mac'][u'srcMac'])

    def indication_parent_change(self, old_parent, new_parent):
        assert old_parent != new_parent

        # allocate the same number of cells to the new parent as it has for the
        # old parent; note that there could be three types of cells:
        # (TX=1,RX=1,SHARED=1), (TX=1), and (RX=1)
        if old_parent is None:
            num_tx_cells = self.NUM_INITIAL_NEGOTIATED_TX_CELLS
            num_rx_cells = self.NUM_INITIAL_NEGOTIATED_RX_CELLS
        else:
            dedicated_cells = self.mote.tsch.get_cells(
                mac_addr         = old_parent,
                slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
            )
            num_tx_cells = len(
                [cell for cell in dedicated_cells if cell.options == [d.CELLOPTION_TX]]
            )
            if num_tx_cells < self.NUM_INITIAL_NEGOTIATED_TX_CELLS:
                num_tx_cells = self.NUM_INITIAL_NEGOTIATED_TX_CELLS
            num_rx_cells = len(
                [cell for cell in dedicated_cells if cell.options == [d.CELLOPTION_RX]]
            )
            if num_rx_cells < self.NUM_INITIAL_NEGOTIATED_RX_CELLS:
                num_rx_cells = self.NUM_INITIAL_NEGOTIATED_RX_CELLS
        if new_parent:
            # reset the retry counter
            # we may better to make sure there is no outstanding
            # transaction with the same peer
            self.retry_count[new_parent] = 0
            self._request_adding_cells(
                neighbor       = new_parent,
                num_tx_cells   = num_tx_cells,
                num_rx_cells   = num_rx_cells
            )

        # clear all the cells allocated for the old parent
        def _callback(event, packet):
            if event == d.SIXP_CALLBACK_EVENT_FAILURE:
                # optimization which is not mentioned in 6P/MSF spec: remove
                # the outstanding transaction because we're deleting all the
                # cells scheduled to the peer now. The outstanding transaction
                # should have the same transaction key as the packet we were
                # trying to send.
                self.mote.sixp.abort_transaction(
                    initiator_mac_addr=packet[u'mac'][u'srcMac'],
                    responder_mac_addr=packet[u'mac'][u'dstMac']
                )
            self._clear_cells(old_parent)

        if old_parent:
            cells = self.mote.tsch.get_cells(
                mac_addr         = old_parent,
                slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
            )
            if len(cells) >= 1:
                self.mote.sixp.send_request(
                    dstMac   = old_parent,
                    command  = d.SIXP_CMD_CLEAR,
                    callback = _callback
                )
            else:
                # do nothing
                pass

    def detect_schedule_inconsistency(self, peerMac):
        # send a CLEAR request to the peer
        self.mote.sixp.send_request(
            dstMac   = peerMac,
            command  = d.SIXP_CMD_CLEAR,
            callback = lambda event, packet: self._clear_cells(peerMac)
        )

    def recv_request(self, packet):
        if   packet[u'app'][u'code'] == d.SIXP_CMD_ADD:
            self._receive_add_request(packet)
        elif packet[u'app'][u'code'] == d.SIXP_CMD_DELETE:
            self._receive_delete_request(packet)
        elif packet[u'app'][u'code'] == d.SIXP_CMD_CLEAR:
            self._receive_clear_request(packet)
        elif packet[u'app'][u'code'] == d.SIXP_CMD_RELOCATE:
            self._receive_relocate_request(packet)
        else:
            # not implemented or not supported
            # ignore this request
            pass

    def clear_to_send_EBs_DATA(self):
        # True if we have a TX cell to the current parent
        slotframe = self.mote.tsch.get_slotframe(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)
        parent_addr = self.mote.rpl.getPreferredParent()
        if (
                (slotframe is None)
                or
                (parent_addr is None)
            ):
            tx_cells = []
        else:
            tx_cells = [
                cell for cell in slotframe.get_cells_by_mac_addr(parent_addr)
                if cell.options == [d.CELLOPTION_TX]
            ]

        if self.mote.dagRoot:
            ret_val = True
        else:
            ret_val = bool(tx_cells)

        return ret_val

    def get_tx_cells(self, mac_addr):
        slotframe = self.mote.tsch.get_slotframe(
            self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
        )
        if slotframe:
            cells = slotframe.get_cells_by_mac_addr(mac_addr)
            negotiated_tx_cells = [cell for cell in cells
                                   if cell.options == [d.CELLOPTION_TX]]
            autonomous_tx_cell = self.get_autonomous_tx_cell(mac_addr)
            if negotiated_tx_cells:
                assert not autonomous_tx_cell
                ret = negotiated_tx_cells
            elif autonomous_tx_cell:
                ret = [autonomous_tx_cell]
            else:
                ret = []
        else:
            ret = []
        return ret

    def get_negotiated_rx_cells(self, mac_addr):
        slotframe = self.mote.tsch.get_slotframe(
            self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
        )
        if slotframe:
            cells = slotframe.get_cells_by_mac_addr(mac_addr)
            ret = [
                cell for cell in cells
                if cell.options == [d.CELLOPTION_RX]
            ]
        else:
            ret = []
        return ret

    def get_autonomous_rx_cell(self):
        slotframe = self.mote.tsch.get_slotframe(
            self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS
        )
        if slotframe:
            cells = slotframe.get_cells_by_mac_addr(None)
            if cells:
                assert len(cells) == 1
                assert cells[0].options == [d.CELLOPTION_RX]
                ret = cells[0]
            else:
                ret = None
        else:
            ret = None
        return ret

    def allocate_autonomous_rx_cell(self):
        mac_addr = self.mote.get_mac_addr()
        slot_offset, channel_offset = self._compute_autonomous_cell(mac_addr)
        self.mote.tsch.addCell(
            slotOffset       = slot_offset,
            channelOffset    = channel_offset,
            neighbor         = None,
            cellOptions      = [
                d.CELLOPTION_RX
            ],
            slotframe_handle = self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS
        )

    def get_autonomous_tx_cell(self, mac_addr):
        slotframe = self.mote.tsch.get_slotframe(
            self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS
        )
        if slotframe:
            cells = slotframe.get_cells_by_mac_addr(mac_addr)
            autonomous_cells = [
                cell for cell in cells
                if (
                        (d.CELLOPTION_TX in cell.options)
                        and
                        (d.CELLOPTION_SHARED in cell.options)
                )
            ]
            if autonomous_cells:
                assert len(autonomous_cells) == 1
                ret = autonomous_cells[0]
            else:
                ret = None
        else:
            ret = None
        return ret

    def allocate_autonomous_tx_cell(self, mac_addr):
        slot_offset, channel_offset = self._compute_autonomous_cell(mac_addr)
        self.mote.tsch.addCell(
            slotOffset       = slot_offset,
            channelOffset    = channel_offset,
            neighbor         = mac_addr,
            cellOptions      = [
                d.CELLOPTION_TX,
                d.CELLOPTION_SHARED
            ],
            slotframe_handle = self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS
        )

    def deallocate_autonomous_tx_cell(self, mac_addr):
        slot_offset, channel_offset = self._compute_autonomous_cell(mac_addr)
        self.mote.tsch.deleteCell(
            slotOffset       = slot_offset,
            channelOffset    = channel_offset,
            neighbor         = mac_addr,
            cellOptions      = [
                d.CELLOPTION_TX,
                d.CELLOPTION_SHARED
            ],
            slotframe_handle = self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS
        )

    # ======================= private ==========================================

    def _reset_cell_counters(self, cell_opt):
        if cell_opt == self.TX_CELL_OPT:
            self.num_tx_cells_elapsed = 0
            self.num_tx_cells_used = 0
        else:
            assert cell_opt == self.RX_CELL_OPT
            self.num_rx_cells_elapsed = 0
            self.num_rx_cells_used = 0

    def _handle_rx_cell_elapsed_event(self, used_by_parent):
        preferred_parent = self.mote.rpl.getPreferredParent()
        self._update_cell_counters(self.RX_CELL_OPT, used_by_parent)
        # adapt number of cells if necessary
        rx_cell_utilization = (
            self.num_rx_cells_used /
            float(self.num_rx_cells_elapsed)
        )
        if d.MSF_MAX_NUMCELLS <= self.num_rx_cells_elapsed:
            if rx_cell_utilization != self.rx_cell_utilization:
                self.log(
                    SimEngine.SimLog.LOG_MSF_RX_CELL_UTILIZATION,
                    {
                        u'_mote_id'    : self.mote.id,
                        u'neighbor'    : preferred_parent,
                        u'value'       : u'{0}% -> {1}%'.format(
                            int(self.rx_cell_utilization * 100),
                            int(rx_cell_utilization * 100)
                        )
                    }
                )
                self.rx_cell_utilization = rx_cell_utilization
            self._adapt_to_traffic(preferred_parent, self.RX_CELL_OPT)
            self._reset_cell_counters(self.RX_CELL_OPT)

    def _update_cell_counters(self, cell_opt, used):
        if cell_opt == self.TX_CELL_OPT:
            self.num_tx_cells_elapsed += 1
            if used:
                self.num_tx_cells_used += 1
        else:
            assert cell_opt == self.RX_CELL_OPT
            self.num_rx_cells_elapsed += 1
            if used:
                self.num_rx_cells_used += 1

    def _adapt_to_traffic(self, neighbor, cell_opt):
        # reset retry counter
        assert neighbor in self.retry_count
        if self.retry_count[neighbor] != -1:
            # we're in the middle of a 6P transaction; try later
            return

        if cell_opt == self.TX_CELL_OPT:
            if d.MSF_LIM_NUMCELLSUSED_HIGH < self.tx_cell_utilization:
                # add one TX cell
                self.retry_count[neighbor] = 0
                self._request_adding_cells(
                    neighbor     = neighbor,
                    num_tx_cells = 1
                )

            elif self.tx_cell_utilization < d.MSF_LIM_NUMCELLSUSED_LOW:
                tx_cells = [cell for cell in self.mote.tsch.get_cells(
                        neighbor,
                        self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                    ) if cell.options == [d.CELLOPTION_TX]]
                # delete one *TX* cell but we need to keep one dedicated
                # cell to our parent at least
                if len(tx_cells) > 1:
                    self.retry_count[neighbor] = 0
                    self._request_deleting_cells(
                        neighbor     = neighbor,
                        num_cells    = 1,
                        cell_options = self.TX_CELL_OPT
                    )
        else:
            assert cell_opt == self.RX_CELL_OPT
            if d.MSF_LIM_NUMCELLSUSED_HIGH < self.rx_cell_utilization:
                self.retry_count[neighbor] = 0
                self._request_adding_cells(
                    neighbor     = neighbor,
                    num_tx_cells = 0,
                    num_rx_cells = 1,
                )

            elif self.rx_cell_utilization < d.MSF_LIM_NUMCELLSUSED_LOW:
                rx_cells = [cell for cell in self.mote.tsch.get_cells(
                        neighbor,
                        self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                    ) if cell.options == [d.CELLOPTION_RX]]
                # delete one *TX* cell but we need to keep one dedicated
                # cell to our parent at least
                if len(rx_cells) > self.NUM_INITIAL_NEGOTIATED_RX_CELLS:
                    self.retry_count[neighbor] = 0
                    self._request_deleting_cells(
                        neighbor     = neighbor,
                        num_cells    = 1,
                        cell_options = self.RX_CELL_OPT
                    )


    def _housekeeping_collision(self):
        """
        Identify cells where schedule collisions occur.
        draft-chang-6tisch-msf-01:
            The key for detecting a schedule collision is that, if a node has
            several cells to the same preferred parent, all cells should exhibit
            the same PDR.  A cell which exhibits a PDR significantly lower than
            the others indicates than there are collisions on that cell.
        :return:
        """

        if self.mote.tsch.get_slotframe(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS) is None:
            return

        # for quick access; get preferred parent
        preferred_parent = self.mote.rpl.getPreferredParent()

        # collect TX cells which has enough numTX
        tx_cell_list = [cell for cell in self.mote.tsch.get_cells(preferred_parent, self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS) if cell.options == [d.CELLOPTION_TX]]
        # pick up TX cells whose NumTx is larger than
        # MSF_MIN_NUM_TX. This is an implementation decision, which is
        # easier to implement than what section 5.3 of
        # draft-ietf-6tisch-msf-03.txt describes as the step-2 of the
        # house-keeping process.
        tx_cell_list = {
            cell.slot_offset: cell for cell in tx_cell_list if (
                d.MSF_MIN_NUM_TX < cell.num_tx
            )
        }
        # collect PDRs of the TX cells
        def pdr(cell):
            assert cell.num_tx > 0
            return cell.num_tx_ack / float(cell.num_tx)
        pdr_list = {
            slotOffset: pdr(cell) for slotOffset, cell in list(tx_cell_list.items())
        }

        if len(pdr_list) > 0:
            # find a cell to relocate using the highest PDR value
            highest_pdr = max(pdr_list.values())
            relocation_cell_list = [
                {
                    'slotOffset'   : slotOffset,
                    'channelOffset': tx_cell_list[slotOffset].channel_offset
                } for slotOffset, pdr in list(pdr_list.items()) if (
                    d.MSF_RELOCATE_PDRTHRES < (highest_pdr - pdr)
                )
            ]
            if (
                    len(relocation_cell_list) > 0
                    and
                    self.retry_count[preferred_parent] == -1
                ):
                # reset retry counter
                self.retry_count[preferred_parent] = 0
                self._request_relocating_cells(
                    neighbor             = preferred_parent,
                    cell_options         = self.TX_CELL_OPT,
                    num_relocating_cells = len(relocation_cell_list),
                    cell_list            = relocation_cell_list
                )
        else:
            # we don't have any TX cell whose PDR is available; do nothing
            pass

        # schedule next housekeeping
        self.engine.scheduleIn(
            delay         = d.MSF_HOUSEKEEPINGCOLLISION_PERIOD,
            cb            = self._housekeeping_collision,
            uniqueTag     = (self.mote.id, u'_housekeeping_collision'),
            intraSlotOrder= d.INTRASLOTORDER_STACKTASKS,
        )

    # cell manipulation helpers
    def _lock_cells(self, cell_list):
        for cell in cell_list:
            self.locked_slots.add(cell[u'slotOffset'])

    def _unlock_cells(self, cell_list):
        for cell in cell_list:
            self.locked_slots.remove(cell[u'slotOffset'])

    def _add_cells(self, neighbor, cell_list, cell_options):
        try:
            for cell in cell_list:
                self.mote.tsch.addCell(
                    slotOffset         = cell[u'slotOffset'],
                    channelOffset      = cell[u'channelOffset'],
                    neighbor           = neighbor,
                    cellOptions        = cell_options,
                    slotframe_handle   = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                )
            if (
                    cell_options == [d.CELLOPTION_TX]
                    and
                    self.get_autonomous_tx_cell(neighbor)
                ):
                self.deallocate_autonomous_tx_cell(neighbor)
        except Exception:
            # We may fail in adding cells since they could be allocated for
            # another peer. We need to have a locking or reservation mechanism
            # to avoid such a situation.
            raise

    def _delete_cells(self, neighbor, cell_list, cell_options):
        for cell in cell_list:
            if self.mote.tsch.get_cell(
                    slot_offset      = cell[u'slotOffset'],
                    channel_offset   = cell[u'channelOffset'],
                    mac_addr         = neighbor,
                    slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
               ) is None:
                # the cell may have been deleted for some reason
                continue
            self.mote.tsch.deleteCell(
                slotOffset       = cell[u'slotOffset'],
                channelOffset    = cell[u'channelOffset'],
                neighbor         = neighbor,
                cellOptions      = cell_options,
                slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
            )

    def _clear_cells(self, neighbor):
        cells = self.mote.tsch.get_cells(
            neighbor,
            self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
        )
        for cell in cells:
            assert neighbor == cell.mac_addr
            assert d.CELLOPTION_SHARED not in cell.options
            self.mote.tsch.deleteCell(
                slotOffset       = cell.slot_offset,
                channelOffset    = cell.channel_offset,
                neighbor         = cell.mac_addr,
                cellOptions      = cell.options,
                slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
            )
        self.mote.sixp.reset_seqnum(neighbor)

    def _relocate_cells(
            self,
            neighbor,
            src_cell_list,
            dst_cell_list,
            cell_options
        ):
        if not dst_cell_list:
            return

        assert len(src_cell_list) == len(dst_cell_list)
        # relocation
        self._add_cells(neighbor, dst_cell_list, cell_options)
        self._delete_cells(neighbor, src_cell_list, cell_options)

    def _get_available_slots(self):
        return list(
            set(self.mote.tsch.get_available_slots(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)) -
            self.locked_slots
        )

    def _create_available_cell_list(self, cell_list_len):
        available_slots = self._get_available_slots()
        # remove slot offset 0 that is reserved for the minimal shared
        # cell
        if 0 in available_slots:
            available_slots.remove(0)

        # remove the slot offset used for the autonomous RX cell
        autonomous_rx_cell = self.get_autonomous_rx_cell()
        assert autonomous_rx_cell
        if autonomous_rx_cell.slot_offset in available_slots:
            available_slots.remove(autonomous_rx_cell.slot_offset)

        if len(available_slots) < cell_list_len:
            # we don't have enough available cells; no cell is selected
            selected_slots = []
        else:
            selected_slots = random.sample(available_slots, cell_list_len)

        cell_list = []
        for slot_offset in selected_slots:
            channel_offset = random.randint(0, self.settings.phy_numChans - 1)
            cell_list.append(
                {
                    'slotOffset'   : slot_offset,
                    'channelOffset': channel_offset
                }
            )
        self._lock_cells(cell_list)
        return cell_list

    def _create_occupied_cell_list(
            self,
            neighbor,
            cell_options,
            cell_list_len
        ):

        occupied_cells = [cell for cell in self.mote.tsch.get_cells(neighbor, self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS) if cell.options == cell_options]

        cell_list = [
            {
                'slotOffset'   : cell.slot_offset,
                'channelOffset': cell.channel_offset
            } for cell in occupied_cells
        ]

        if cell_list_len <= len(occupied_cells):
            cell_list = random.sample(cell_list, cell_list_len)

        return cell_list

    def _are_cells_allocated(
            self,
            peerMac,
            cell_list,
            cell_options
        ):

        # collect allocated cells
        assert cell_options in [self.TX_CELL_OPT, self.RX_CELL_OPT]
        allocated_cells = [cell for cell in self.mote.tsch.get_cells(peerMac, self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS) if cell.options == cell_options]

        # test all the cells in the cell list against the allocated cells
        ret_val = True
        for cell in cell_list:
            slotOffset    = cell[u'slotOffset']
            channelOffset = cell[u'channelOffset']
            cell = self.mote.tsch.get_cell(
                slot_offset      = slotOffset,
                channel_offset   = channelOffset,
                mac_addr         = peerMac,
                slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
            )

            if cell is None:
                ret_val = False
                break

        return ret_val

    # ADD command related stuff
    def _request_adding_cells(
            self,
            neighbor,
            num_tx_cells   = 0,
            num_rx_cells   = 0
        ):

        # determine num_cells and cell_options; update num_{tx,rx}_cells
        if num_tx_cells > 0:
            cell_options = self.TX_CELL_OPT
            if num_tx_cells < self.DEFAULT_CELL_LIST_LEN:
                num_cells    = num_tx_cells
                num_tx_cells = 0
            else:
                num_cells    = self.DEFAULT_CELL_LIST_LEN
                num_tx_cells = num_tx_cells - self.DEFAULT_CELL_LIST_LEN
        elif num_rx_cells > 0:
            cell_options = self.RX_CELL_OPT
            num_cells    = num_rx_cells
            if num_rx_cells < self.DEFAULT_CELL_LIST_LEN:
                num_cells    = num_rx_cells
                num_rx_cells = 0
            else:
                num_cells    = self.DEFAULT_CELL_LIST_LEN
                num_rx_cells = num_rx_cells - self.DEFAULT_CELL_LIST_LEN
        else:
            # nothing to add
            self.retry_count[neighbor] = -1
            return

        # prepare cell_list
        cell_list = self._create_available_cell_list(self.DEFAULT_CELL_LIST_LEN)

        if len(cell_list) == 0:
            # we don't have available cells right now
            self.log(
                SimEngine.SimLog.LOG_MSF_ERROR_SCHEDULE_FULL,
                {
                    '_mote_id'    : self.mote.id
                }
            )
            self.retry_count[neighbor] = -1
            return

        # prepare _callback which is passed to SixP.send_request()
        callback = self._create_add_request_callback(
            neighbor,
            num_cells,
            cell_options,
            cell_list,
            num_tx_cells,
            num_rx_cells
        )

        # send a request
        self.mote.sixp.send_request(
            dstMac      = neighbor,
            command     = d.SIXP_CMD_ADD,
            cellOptions = cell_options,
            numCells    = num_cells,
            cellList    = cell_list,
            callback    = callback
        )

    def _receive_add_request(self, request):

        # for quick access
        proposed_cells = request[u'app'][u'cellList']
        peerMac         = request[u'mac'][u'srcMac']

        # find available cells in the received CellList
        slots_in_cell_list = set(
            [c[u'slotOffset'] for c in proposed_cells]
        )
        available_slots  = list(
            slots_in_cell_list.intersection(
                set(self._get_available_slots())
            )
        )

        # prepare cell_list
        candidate_cells = [
            c for c in proposed_cells if c[u'slotOffset'] in available_slots
        ]
        if len(candidate_cells) < request[u'app'][u'numCells']:
            cell_list = candidate_cells
        else:
            cell_list = random.sample(
                candidate_cells,
                request[u'app'][u'numCells']
            )

        # prepare callback
        if len(available_slots) > 0:
            code = d.SIXP_RC_SUCCESS

            self._lock_cells(candidate_cells)
            def callback(event, packet):
                if event == d.SIXP_CALLBACK_EVENT_MAC_ACK_RECEPTION:
                    # prepare cell options for this responder
                    if request[u'app'][u'cellOptions'] == self.TX_CELL_OPT:
                        # invert direction
                        cell_options = self.RX_CELL_OPT
                    elif request[u'app'][u'cellOptions'] == self.RX_CELL_OPT:
                        # invert direction
                        cell_options = self.TX_CELL_OPT
                    else:
                        # Unsupported cell options for MSF
                        raise Exception()

                    self._add_cells(
                        neighbor     = peerMac,
                        cell_list    = cell_list,
                        cell_options = cell_options
                )
                self._unlock_cells(candidate_cells)
        else:
            code      = d.SIXP_RC_ERR
            cell_list = None
            callback  = None

        # send a response
        self.mote.sixp.send_response(
            dstMac      = peerMac,
            return_code = code,
            cellList    = cell_list,
            callback    = callback
        )

    def _create_add_request_callback(
            self,
            neighbor,
            num_cells,
            cell_options,
            cell_list,
            num_tx_cells,
            num_rx_cells
        ):
        def callback(event, packet):
            if event == d.SIXP_CALLBACK_EVENT_PACKET_RECEPTION:
                assert packet[u'app'][u'msgType'] == d.SIXP_MSG_TYPE_RESPONSE
                if packet[u'app'][u'code'] == d.SIXP_RC_SUCCESS:
                    # add cells on success of the transaction
                    self._add_cells(
                        neighbor     = neighbor,
                        cell_list    = packet[u'app'][u'cellList'],
                        cell_options = cell_options
                    )
                    # The received CellList could be smaller than the requested
                    # NumCells; adjust num_{tx,rx}_cells
                    _num_tx_cells   = num_tx_cells
                    _num_rx_cells   = num_rx_cells
                    remaining_cells = num_cells - len(packet[u'app'][u'cellList'])
                    if remaining_cells > 0:
                        if cell_options == self.TX_CELL_OPT:
                            _num_tx_cells -= remaining_cells
                        elif cell_options == self.RX_CELL_OPT:
                            _num_rx_cells -= remaining_cells
                        else:
                            # never comes here
                            raise Exception()

                    # start another transaction
                    self.retry_count[neighbor] = 0
                    self._request_adding_cells(
                        neighbor       = neighbor,
                        num_tx_cells   = _num_tx_cells,
                        num_rx_cells   = _num_rx_cells
                    )
                    #print(str(self.mote.tsch.get_available_slots(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)))
                    #print(str(cell_list))
                else:
                    # TODO: request doesn't succeed; how should we do?
                    self.retry_count[neighbor] = -1

            elif event == d.SIXP_CALLBACK_EVENT_TIMEOUT:
                if self.retry_count[neighbor] == self.MAX_RETRY:
                    # give up this neighbor
                    if neighbor == self.mote.rpl.getPreferredParent():
                        self.mote.rpl.of.poison_rpl_parent(neighbor)
                    self.retry_count[neighbor] = -1 # done
                else:
                    # retry
                    self.retry_count[neighbor] += 1
                    if cell_options == self.TX_CELL_OPT:
                        _num_tx_cells = num_cells + num_tx_cells
                        _num_rx_cells = num_rx_cells
                    else:
                        _num_tx_cells = num_tx_cells
                        _num_rx_cells = num_cells + num_rx_cells
                    self._request_adding_cells(
                        neighbor       = neighbor,
                        num_tx_cells   = _num_tx_cells,
                        num_rx_cells   = _num_rx_cells
                    )
            else:
                # ignore other events
                pass

            # unlock the slots used in this transaction
            self._unlock_cells(cell_list)

        return callback

    # DELETE command related stuff
    def _request_deleting_cells(
            self,
            neighbor,
            num_cells,
            cell_options
        ):

        # prepare cell_list to send
        cell_list = self._create_occupied_cell_list(
            neighbor      = neighbor,
            cell_options  = cell_options,
            cell_list_len = self.DEFAULT_CELL_LIST_LEN
        )
        assert len(cell_list) > 0

        # prepare callback
        callback = self._create_delete_request_callback(
            neighbor,
            num_cells,
            cell_options
        )

        # send a DELETE request
        self.mote.sixp.send_request(
            dstMac      = neighbor,
            command     = d.SIXP_CMD_DELETE,
            cellOptions = cell_options,
            numCells    = num_cells,
            cellList    = cell_list,
            callback    = callback
        )

    def _receive_delete_request(self, request):
        # for quick access
        num_cells           = request[u'app'][u'numCells']
        cell_options        = request[u'app'][u'cellOptions']
        candidate_cell_list = request[u'app'][u'cellList']
        peerMac             = request[u'mac'][u'srcMac']

        # confirm all the cells in the cell list are allocated for the peer
        # with the specified cell options
        #
        # invert the direction in cell_options
        assert cell_options in [self.TX_CELL_OPT, self.RX_CELL_OPT]
        if   cell_options == self.TX_CELL_OPT:
            our_cell_options = self.RX_CELL_OPT
        elif cell_options == self.RX_CELL_OPT:
            our_cell_options   = self.TX_CELL_OPT

        if (
                (
                    self._are_cells_allocated(
                        peerMac      = peerMac,
                        cell_list    = candidate_cell_list,
                        cell_options = our_cell_options
                    ) is True
                )
                and
                (num_cells <= len(candidate_cell_list))
            ):
            code = d.SIXP_RC_SUCCESS
            cell_list = random.sample(candidate_cell_list, num_cells)

            def callback(event, packet):
                if event == d.SIXP_CALLBACK_EVENT_MAC_ACK_RECEPTION:
                    self._delete_cells(
                        neighbor     = peerMac,
                        cell_list    = cell_list,
                        cell_options = our_cell_options
                )
        else:
            code      = d.SIXP_RC_ERR
            cell_list = None
            callback  = None

        # send the response
        self.mote.sixp.send_response(
            dstMac      = peerMac,
            return_code = code,
            cellList    = cell_list,
            callback    = callback
        )

    def _create_delete_request_callback(
            self,
            neighbor,
            num_cells,
            cell_options
        ):
        def callback(event, packet):
            if (
                    (event == d.SIXP_CALLBACK_EVENT_PACKET_RECEPTION)
                    and
                    (packet[u'app'][u'msgType'] == d.SIXP_MSG_TYPE_RESPONSE)
                ):
                self.retry_count[neighbor] = -1
                if packet[u'app'][u'code'] == d.SIXP_RC_SUCCESS:
                    self._delete_cells(
                        neighbor     = neighbor,
                        cell_list    = packet[u'app'][u'cellList'],
                        cell_options = cell_options
                    )
                else:
                    # TODO: request doesn't succeed; how should we do?
                    pass
            elif event == d.SIXP_CALLBACK_EVENT_TIMEOUT:
                if self.retry_count[neighbor] == self.MAX_RETRY:
                    # give it up
                    self.retry_count[neighbor] = -1
                    if neighbor == self.mote.rpl.getPreferredParent():
                        self.mote.rpl.of.poison_rpl_parent(neighbor)
                else:
                    # retry
                    self.retry_count[neighbor] += 1
                    self._request_deleting_cells(
                        neighbor,
                        num_cells,
                        cell_options
                    )
            else:
                # ignore other events
                pass

        return callback

    # RELOCATE command related stuff
    def _request_relocating_cells(
            self,
            neighbor,
            cell_options,
            num_relocating_cells,
            cell_list
        ):

        # determine num_cells and relocation_cell_list;
        # update num_relocating_cells and cell_list
        if self.DEFAULT_CELL_LIST_LEN < num_relocating_cells:
            num_cells             = self.DEFAULT_CELL_LIST_LEN
            relocation_cell_list  = cell_list[:self.DEFAULT_CELL_LIST_LEN]
            num_relocating_cells -= self.DEFAULT_CELL_LIST_LEN
            cell_list             = cell_list[self.DEFAULT_CELL_LIST_LEN:]
        else:
            num_cells             = num_relocating_cells
            relocation_cell_list  = cell_list
            num_relocating_cells  = 0
            cell_list             = []

        # we don't have any cell to relocate; done
        if len(relocation_cell_list) == 0:
            self.retry_count[neighbor] = -1
            return

        # prepare candidate_cell_list
        candidate_cell_list = self._create_available_cell_list(
            self.DEFAULT_CELL_LIST_LEN
        )

        if len(candidate_cell_list) == 0:
            # no available cell to move the cells to
            self.log(
                SimEngine.SimLog.LOG_MSF_ERROR_SCHEDULE_FULL,
                {
                    '_mote_id'    : self.mote.id
                }
            )
            self.retry_count[neighbor] = -1
            return

        # prepare callback
        def callback(event, packet):
            if event == d.SIXP_CALLBACK_EVENT_PACKET_RECEPTION:
                assert packet[u'app'][u'msgType'] == d.SIXP_MSG_TYPE_RESPONSE
                if packet[u'app'][u'code'] == d.SIXP_RC_SUCCESS:
                    # perform relocations
                    num_relocations = len(packet[u'app'][u'cellList'])
                    self._relocate_cells(
                        neighbor      = neighbor,
                        src_cell_list = relocation_cell_list[:num_cells],
                        dst_cell_list = packet[u'app'][u'cellList'],
                        cell_options  = cell_options
                    )

                    # adjust num_relocating_cells and cell_list
                    _num_relocating_cells = (
                        num_relocating_cells + num_cells - num_relocations
                    )
                    _cell_list = (
                        cell_list + relocation_cell_list[num_relocations:]
                    )

                    # start another transaction
                    self.retry_count[neighbor] = 0
                    self._request_relocating_cells(
                        neighbor             = neighbor,
                        cell_options         = cell_options,
                        num_relocating_cells = _num_relocating_cells,
                        cell_list            = _cell_list
                    )
            elif event == d.SIXP_CALLBACK_EVENT_TIMEOUT:
                if self.retry_count[neighbor] == self.MAX_RETRY:
                    # give up this neighbor
                    if neighbor == self.mote.rpl.getPreferredParent():
                        self.mote.rpl.of.poison_rpl_parent(neighbor)
                    self.retry_count[neighbor] = -1 # done
                else:
                    # retry
                    self.retry_count[neighbor] += 1
                    self._request_relocating_cells(
                        neighbor,
                        cell_options,
                        num_relocating_cells,
                        cell_list
                    )

            # unlock the slots used in this transaction
            self._unlock_cells(candidate_cell_list)

        # send a request
        self.mote.sixp.send_request(
            dstMac             = neighbor,
            command            = d.SIXP_CMD_RELOCATE,
            cellOptions        = cell_options,
            numCells           = num_cells,
            relocationCellList = relocation_cell_list,
            candidateCellList  = candidate_cell_list,
            callback           = callback
        )

    def _receive_relocate_request(self, request):
        # for quick access
        num_cells        = request[u'app'][u'numCells']
        cell_options     = request[u'app'][u'cellOptions']
        relocating_cells = request[u'app'][u'relocationCellList']
        candidate_cells  = request[u'app'][u'candidateCellList']
        peerMac          = request[u'mac'][u'srcMac']

        # confirm all the cells in the cell list are allocated for the peer
        # with the specified cell options
        #
        # invert the direction in cell_options
        assert cell_options in [self.TX_CELL_OPT, self.RX_CELL_OPT]
        if   cell_options == self.TX_CELL_OPT:
            our_cell_options = self.RX_CELL_OPT
        elif cell_options == self.RX_CELL_OPT:
            our_cell_options   = self.TX_CELL_OPT

        if (
                (
                    self._are_cells_allocated(
                        peerMac      = peerMac,
                        cell_list    = relocating_cells,
                        cell_options = our_cell_options
                    ) is True
                )
                and
                (num_cells <= len(candidate_cells))
            ):
            # find available cells in the received candidate cell list
            slots_in_slotframe = set(range(0, self.settings.tsch_slotframeLength))
            slots_in_use       = set(
                self.mote.tsch.get_busy_slots(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)
            )
            candidate_slots    = set(
                [c[u'slotOffset'] for c in candidate_cells]
            )
            available_slots    = list(
                candidate_slots.intersection(
                    set(self._get_available_slots())
                )
            )

            code = d.SIXP_RC_SUCCESS
            cell_list = []
            if available_slots:
                # prepare response
                selected_slots = random.sample(available_slots, num_cells)
                for cell in candidate_cells:
                    if cell[u'slotOffset'] in selected_slots:
                        cell_list.append(cell)

                self._lock_cells(cell_list)
            else:
                # we will return an empty cell list with RC_SUCCESS
                pass

            # prepare callback
            def callback(event, packet):
                if event == d.SIXP_CALLBACK_EVENT_MAC_ACK_RECEPTION:
                    num_relocations = len(cell_list)
                    self._relocate_cells(
                        neighbor      = peerMac,
                        src_cell_list = relocating_cells[:num_relocations],
                        dst_cell_list = cell_list,
                        cell_options  = our_cell_options
                    )
                self._unlock_cells(cell_list)

        else:
            code      = d.SIXP_RC_ERR
            cell_list = None
            callback  = None

        # send a response
        self.mote.sixp.send_response(
            dstMac      = peerMac,
            return_code = code,
            cellList    = cell_list,
            callback    = callback
        )


    # CLEAR command related stuff
    def _receive_clear_request(self, request):

        peerMac = request[u'mac'][u'srcMac']

        def callback(event, packet):
            # remove all the cells no matter what happens
            self._clear_cells(peerMac)

        # create CLEAR response
        self.mote.sixp.send_response(
            dstMac      = peerMac,
            return_code = d.SIXP_RC_SUCCESS,
            callback    = callback
        )

    # autonomous cell
    def _compute_autonomous_cell(self, mac_addr):
        slotframe = self.mote.tsch.get_slotframe(
            self.SLOTFRAME_HANDLE_AUTONOMOUS_CELLS
        )
        hash_value = self._sax(mac_addr)

        slot_offset = int(1 + (hash_value % (slotframe.length - 1)))
        channel_offset = int(hash_value % self.settings.phy_numChans)

        return (slot_offset, channel_offset)

    # SAX
    def _sax(self, mac_addr):
        # XXX: a concrete definition of this hash function is needed to be
        # provided by the draft

        LEFT_SHIFT_NUM = 5
        RIGHT_SHIFT_NUM = 2

        # assuming v (seed) is 0
        hash_value = 0
        for word in netaddr.EUI(mac_addr).words:
            for byte in divmod(word, 0x100):
                left_shifted = (hash_value << LEFT_SHIFT_NUM)
                right_shifted = (hash_value >> RIGHT_SHIFT_NUM)
                hash_value ^= left_shifted + right_shifted + byte

        # assuming T (table size) is 16-bit
        return hash_value & 0xFFFF
################################################################################################################################
    
    def DCAB(self, old_preferred=None, new_preferred=None):
        
        self.SlotRange = 101
        self.ChannelRange = 16
        
        motes = self.engine.motes
        root = motes[0]
        tree = root.rpl.tree
        
        treeDict = self.create_tree_Dict(tree)
        schedulingMx = root.sf.schedulingMatrix
        otherParents = root.rpl.otherParents
        packetDroppedPL = root.sf.packetDroppedPerLink
          
        for slot in schedulingMx:
            for channel in schedulingMx[slot]:
                sender_id = schedulingMx[slot][channel][0]
                sender_mote = self.engine.get_mote_by_id(sender_id)
                receiver_id = schedulingMx[slot][channel][1]
                receiver_mote = self.engine.get_mote_by_id(receiver_id)
                sender_mote.sf.retry_count[receiver_mote.get_mac_addr()] = -1
                sender_mote.tsch.deleteCell(
                    slotOffset       = slot,
                    channelOffset    = channel,
                    neighbor         = receiver_mote.get_mac_addr(),
                    cellOptions      = self.TX_CELL_OPT,
                    slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                )
                receiver_mote.tsch.deleteCell(
                    slotOffset       = slot,
                    channelOffset    = channel,
                    neighbor         = sender_mote.get_mac_addr(),
                    cellOptions      = self.RX_CELL_OPT,
                    slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                )
                if len(schedulingMx[slot][channel]) > 2:
                    for parent_id in schedulingMx[slot][channel][2:]:
                        receiver_mote = self.engine.get_mote_by_id(parent_id)
                        receiver_mote.tsch.deleteCell(
                            slotOffset       = slot,
                            channelOffset    = channel,
                            neighbor         = sender_mote.get_mac_addr(),
                            cellOptions      = self.RX_CELL_OPT,
                            slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                        )
                        print("delete(from SF) cells: [",slot,channel,"] for mote ",sender_id,"->",receiver_mote.id)
            
        schedulingMx.clear()
        root.sf.schedulingMatrix = self.create_scheduling_matrix(treeDict,schedulingMx,otherParents)
        schedulingMatrix = root.sf.schedulingMatrix
        
        # print("schedulingMatrix before tentative cells:")        
        # print(schedulingMatrix) 
        
        #schedulingMatrix = self.handle_tentative_cells(treeDict,packetDroppedPL,schedulingMatrix)
        packetDroppedPL.clear()    
        
        # print("schedulingMatrix after tentative cells:")        
        # print(schedulingMatrix)
        self.print_scheduling_matrix(schedulingMatrix)
        for slot in schedulingMatrix:
            for channel in schedulingMatrix[slot]:
                sender_id = schedulingMatrix[slot][channel][0]
                receiver_id = schedulingMatrix[slot][channel][1]
                sender_mote = self.engine.get_mote_by_id(sender_id)
                receiver_mote = self.engine.get_mote_by_id(receiver_id)
                if new_preferred:
                    sender_mote.sf.retry_count[receiver_mote.get_mac_addr()] = 0
                sender_mote.tsch.addCell(
                    slotOffset         = slot,
                    channelOffset      = channel,
                    neighbor           = receiver_mote.get_mac_addr(),
                    cellOptions        = self.TX_CELL_OPT,
                    slotframe_handle   = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                )
                receiver_mote.tsch.addCell(
                    slotOffset         = slot,
                    channelOffset      = channel,
                    neighbor           = sender_mote.get_mac_addr(),
                    cellOptions        = self.RX_CELL_OPT,
                    slotframe_handle   = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                )
                if len(schedulingMatrix[slot][channel]) > 2:
                    for parent_id in schedulingMatrix[slot][channel][2:]:
                        receiver_mote = self.engine.get_mote_by_id(parent_id)
                        receiver_mote.tsch.addCell(
                            slotOffset         = slot,
                            channelOffset      = channel,
                            neighbor           = sender_mote.get_mac_addr(),
                            cellOptions        = self.RX_CELL_OPT,
                            slotframe_handle   = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                        )
                        print("add cells: [",slot,channel,"] for mote ",sender_id,"->",receiver_mote.id)
    
    def handle_tentative_cells(self,treeDict,packetDroppedPL,schedulingMatrix):
        for routes in treeDict:
            for route in treeDict[routes]:
                for link in route:
                    slot = 1
                    child_mote_id = link[0]
                    parent_mote_id = link[1]
                    if child_mote_id in packetDroppedPL:
                        if parent_mote_id in packetDroppedPL[child_mote_id]:
                            tentative_cells = packetDroppedPL[child_mote_id][parent_mote_id]
                            sc = False #just for initializing
                            child_mote = self.engine.get_mote_by_id(child_mote_id)
                            parent_mote = self.engine.get_mote_by_id(parent_mote_id)
                            
                            child_available_slots = child_mote.tsch.get_available_slots(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)
                            parent_available_slots = parent_mote.tsch.get_available_slots(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)
                            
                            #print(child_mote_id,child_available_slots)
                            #print(parent_mote_id,parent_available_slots)

                            if 0 in child_available_slots:
                                child_available_slots.remove(0)
                            if 0 in parent_available_slots:
                                parent_available_slots.remove(0)
                            
                            child_autonomous_rx_cell = child_mote.sf.get_autonomous_rx_cell()
                            parent_autonomous_rx_cell = parent_mote.sf.get_autonomous_rx_cell()
                            
                            assert child_autonomous_rx_cell
                            assert parent_autonomous_rx_cell
                            
                            # print(child_mote_id,child_autonomous_rx_cell)
                            # print(parent_mote_id,parent_autonomous_rx_cell)
                            if child_autonomous_rx_cell.slot_offset in child_available_slots:
                                child_available_slots.remove(child_autonomous_rx_cell.slot_offset)
                            if parent_autonomous_rx_cell.slot_offset in parent_available_slots:
                                parent_available_slots.remove(parent_autonomous_rx_cell.slot_offset)
                            
                            available_slots = list(set(child_available_slots)& set(parent_available_slots))
                            #print(available_slots)
                            for tentative in range(tentative_cells):
                                while slot < 101:
                                    if slot in available_slots:
                                        sc = True
                                        if slot in schedulingMatrix:
                                            ch = 0
                                            for channel in schedulingMatrix[slot]:
                                                if link[0] in schedulingMatrix[slot][channel] or link[1] in schedulingMatrix[slot][channel]:
                                                    sc = False
                                                    slot += 1
                                                    break
                                                else:
                                                    ch = channel+1
                                            if sc:
                                                if ch == self.ChannelRange:
                                                    sc = False
                                                    slot += 1
                                                    continue
                                                else:
                                                    assert ch < self.ChannelRange
                                                    sc = True
                                                    schedulingMatrix[slot][ch] = [link[0],link[1]]
                                                    slot += 1
                                                    break
                                        else:
                                            assert sc == True
                                            sc = True
                                            schedulingMatrix[slot]={}
                                            schedulingMatrix[slot][0] = [link[0],link[1]]
                                            slot += 1
                                            break
                                    else:
                                        slot += 1
                                if not sc:
                                    already_set = False
                                    for slot in schedulingMatrix:
                                        if link in schedulingMatrix[slot].values():
                                            already_set = True
                                            print("there is no slot-channel for link ",str(link)," but this link already sat in scheduling matrix")
                                            break
                                    if not already_set:
                                        print("there is no slot-channel for link ",str(link))
                                        raise SystemExit
        return schedulingMatrix        

    def create_tree_Dict(self,tree):
        # {
            # 1:0,
            # 2:0,
            # 3:0,
            # 4:0,
            # 6:2,
            # 7:1,
            # 8:1,
            # 9:2,
            # 5:8
        # }
        list_dict = {}
        for child in tree:
            parent = tree[child]
            temp_counter = 1
            temp_list = [[child,parent]]
            while parent in tree:
                temp_counter += 1
                temp_list.append([parent,tree[parent]])
                parent = tree[parent]
            #temp_list.reverse()
            if "len{0}".format(temp_counter) in list_dict:
                list_dict["len{0}".format(temp_counter)].append(temp_list)
            else:
                list_dict["len{0}".format(temp_counter)] = []
                list_dict["len{0}".format(temp_counter)].append(temp_list)
        # print(list_dict)
        # {
            # 'len3': [[[5, 8], [8, 1], [1, 0]]],
            # 'len2': [[[6, 2], [2, 0]], [[7, 1], [1, 0]], [[8, 1], [1, 0]], [[9, 2], [2, 0]]],
            # 'len1': [[[1, 0]], [[2, 0]], [[3, 0]], [[4, 0]]]
        # }
        return dict(reversed(sorted(list_dict.items())))
    
    def create_scheduling_matrix(self,treeDict,schedulingMatrix,otherParents):
        # {
            # 1:{
                    # 0:[5,8],
                    # 1:[6,2],
                    # 2:[7,1],
                    # 3:[3,0]
              # },
            # 2:{
                    # 0:[8,1],
                    # 1:[2,0]
              # },
            # 3:{
                    # 0:[1,0],
                    # 1:[9,2]
              # },
            # 4:{
                    # 0:[1,0],
              # },
            # 5:{
                    # 0:[8,1],
                    # 1:[2,0]
              # },
            # 6:{
                    # 0:[1,0],
              # },
            # 7:{
                    # 0:[1,0],
              # },
            # 8:{
                    # 0:[2,0],
              # },
            # 9:{
                    # 0:[4,0],
              # }
        # }
        # {
            # 'len3': [[[5, 8], [8, 1], [1, 0]]],
            # 'len2': [[[6, 2], [2, 0]], [[7, 1], [1, 0]], [[8, 1], [1, 0]], [[9, 2], [2, 0]]],
            # 'len1': [[[1, 0]], [[2, 0]], [[3, 0]], [[4, 0]]]
        # }
        for routes in treeDict:
            for route in treeDict[routes]:
                slot = 1
                for link in route:
                    sc = False #just for initializing
                    child_mote_id = link[0]
                    parent_mote_id = link[1]
                    print(child_mote_id,parent_mote_id)
                    child_mote = self.engine.get_mote_by_id(child_mote_id)
                    parent_mote = self.engine.get_mote_by_id(parent_mote_id)
                    
                    child_available_slots = child_mote.tsch.get_available_slots(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)
                    parent_available_slots = parent_mote.tsch.get_available_slots(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)
                    
                    #print(child_mote_id,child_available_slots)
                    #print(parent_mote_id,parent_available_slots)

                    if 0 in child_available_slots:
                        child_available_slots.remove(0)
                    if 0 in parent_available_slots:
                        parent_available_slots.remove(0)
                    
                    child_autonomous_rx_cell = child_mote.sf.get_autonomous_rx_cell()
                    parent_autonomous_rx_cell = parent_mote.sf.get_autonomous_rx_cell()
                    
                    assert child_autonomous_rx_cell
                    assert parent_autonomous_rx_cell
                    
                    # print(child_mote_id,child_autonomous_rx_cell)
                    # print(parent_mote_id,parent_autonomous_rx_cell)
                    if child_autonomous_rx_cell.slot_offset in child_available_slots:
                        child_available_slots.remove(child_autonomous_rx_cell.slot_offset)
                    if parent_autonomous_rx_cell.slot_offset in parent_available_slots:
                        parent_available_slots.remove(parent_autonomous_rx_cell.slot_offset)
                    
                    available_slots = list(set(child_available_slots)& set(parent_available_slots))
                    
                    if child_mote_id in otherParents:
                        print("available_slots before handle another parent")
                        print(available_slots)
                        for another_parent_mote_id in otherParents[child_mote_id]:
                            another_parent_mote = self.engine.get_mote_by_id(another_parent_mote_id)
                            another_parent_available_slots = another_parent_mote.tsch.get_available_slots(self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS)
                            if 0 in another_parent_available_slots:
                                another_parent_available_slots.remove(0)
                            another_parent_autonomous_rx_cell = another_parent_mote.sf.get_autonomous_rx_cell()
                            if another_parent_autonomous_rx_cell.slot_offset in another_parent_available_slots:
                                another_parent_available_slots.remove(another_parent_autonomous_rx_cell.slot_offset)
                            available_slots = list(set(available_slots)& set(another_parent_available_slots))
                        print("available_slots after handle another parent")
                        print(available_slots)
                    
                    cell_list = [link[0],link[1]]
                    if link[0] in otherParents:
                        for parent_id in otherParents[link[0]]:
                            cell_list.append(parent_id)
                    while slot < 101:
                        if slot in available_slots:
                            sc = True
                            if slot in schedulingMatrix:
                                ch = 0
                                for channel in schedulingMatrix[slot]:
                                    for mote_id in cell_list:
                                        if mote_id in schedulingMatrix[slot][channel]:
                                            sc = False
                                            slot += 1
                                            break
                                    if sc:
                                        ch = channel+1
                                    else:
                                        break
                                if sc:
                                    if ch == self.ChannelRange:
                                        sc = False
                                        slot += 1
                                        continue
                                    else:
                                        assert ch < self.ChannelRange
                                        sc = True
                                        schedulingMatrix[slot][ch] = cell_list
                                        break
                            else:
                                assert sc == True
                                sc = True
                                schedulingMatrix[slot]={}
                                schedulingMatrix[slot][0] = cell_list
                                break
                        else:
                            slot += 1
                    if not sc:
                        already_set = False
                        for slot in schedulingMatrix:
                            if cell_list in schedulingMatrix[slot].values():
                                already_set = True
                                print("there is no slot-channel for link ",str(link)," but this link already sat in scheduling matrix")
                                break
                        if not already_set:
                            print("there is no slot-channel for link ",str(link))
                            raise SystemExit
        return schedulingMatrix
    
    def delete_old_parents(self,joined_mote,otherParents,tree):
        # called by root in rpl
        schedulingMx = self.schedulingMatrix
        for slot in schedulingMx:
            for channel in schedulingMx[slot]:
                sender_id = schedulingMx[slot][channel][0]
                receiver_id = schedulingMx[slot][channel][1]
                if sender_id == joined_mote.id:
                    sender_mote = self.engine.get_mote_by_id(sender_id)
                    receiver_mote = self.engine.get_mote_by_id(receiver_id)
                    sender_mote.sf.retry_count[receiver_mote.get_mac_addr()] = -1
                    sender_mote.tsch.deleteCell(
                        slotOffset       = slot,
                        channelOffset    = channel,
                        neighbor         = receiver_mote.get_mac_addr(),
                        cellOptions      = self.TX_CELL_OPT,
                        slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                    )
                    receiver_mote.tsch.deleteCell(
                        slotOffset       = slot,
                        channelOffset    = channel,
                        neighbor         = sender_mote.get_mac_addr(),
                        cellOptions      = self.RX_CELL_OPT,
                        slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                    )
                    print("delete(from RPL) cells: [",slot,channel,"] for mote ",sender_id,"->",receiver_mote.id)
                    if len(schedulingMx[slot][channel]) > 2:
                        for parent_id in schedulingMx[slot][channel][2:]:
                            receiver_mote = self.engine.get_mote_by_id(parent_id)
                            receiver_mote.tsch.deleteCell(
                                slotOffset       = slot,
                                channelOffset    = channel,
                                neighbor         = sender_mote.get_mac_addr(),
                                cellOptions      = self.RX_CELL_OPT,
                                slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                            )
                            print("delete(from RPL) cells: [",slot,channel,"] for mote ",sender_id,"->",receiver_mote.id)
                    del schedulingMx[slot][channel]
                    break
        if joined_mote.id in otherParents:
            del otherParents[joined_mote.id]
        if joined_mote.id in tree:
            del tree[joined_mote.id]
        
    def handle_desynced_mote(self,desynced_mote):
        motes = self.engine.motes
        root = motes[0]
        tree = root.rpl.tree
        schedulingMatrix = root.sf.schedulingMatrix
        otherParents = root.rpl.otherParents
        print("tree before desynced")
        print(tree)
        print("otherParents before desynced")
        print(otherParents)
        for slot in schedulingMatrix:
            for channel in schedulingMatrix[slot]:
                sender_id = schedulingMatrix[slot][channel][0]
                receiver_id = schedulingMatrix[slot][channel][1]
                if sender_id == desynced_mote.id:
                    sender_mote = self.engine.get_mote_by_id(sender_id)
                    receiver_mote = self.engine.get_mote_by_id(receiver_id)
                    sender_mote.tsch.deleteCell(
                        slotOffset       = slot,
                        channelOffset    = channel,
                        neighbor         = receiver_mote.get_mac_addr(),
                        cellOptions      = self.TX_CELL_OPT,
                        slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                    )
                    receiver_mote.tsch.deleteCell(
                        slotOffset       = slot,
                        channelOffset    = channel,
                        neighbor         = sender_mote.get_mac_addr(),
                        cellOptions      = self.RX_CELL_OPT,
                        slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                    )
                    print("delete(from desynced-primary parent) cells: [",slot,channel,"] for mote ",sender_id,"->",receiver_mote.id)
                    if len(schedulingMatrix[slot][channel]) > 2:
                        for parent_id in schedulingMatrix[slot][channel][2:]:
                            receiver_mote = self.engine.get_mote_by_id(parent_id)
                            receiver_mote.tsch.deleteCell(
                                slotOffset       = slot,
                                channelOffset    = channel,
                                neighbor         = sender_mote.get_mac_addr(),
                                cellOptions      = self.RX_CELL_OPT,
                                slotframe_handle = self.SLOTFRAME_HANDLE_NEGOTIATED_CELLS
                            )
                            print("delete(from desynced-primary parent) cells: [",slot,channel,"] for mote ",sender_id,"->",receiver_mote.id)
                    del schedulingMatrix[slot][channel]
                    break
        if desynced_mote.id in otherParents:
            del otherParents[desynced_mote.id]
        if desynced_mote.id in tree:
            del tree[desynced_mote.id]
        print("tree after desynced")
        print(tree)
        print("otherParents after desynced")
        print(otherParents)
    
    def print_scheduling_matrix(self,sm):
        print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        print("         slot 1    slot 2    slot 3    slot 4    slot 5    slot 6    slot 7    slot 8    slot 9    slot 10   slot 11   slot 12   slot 13   slot 14   slot 15   slot 16   slot 17   slot 18   slot 19")
        print("       ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        for ch in range(8):
            print "ch "+str(ch),
            for slot in range(1,20):
                if slot in sm.keys():
                    if ch in sm[slot]:
                        print " | "+str(sm[slot][ch]),
                    else:
                        print " |      ",
                else:
                    print " |       ",
            print "\n"
        print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")