SUBROUTINE PAPIF_library_init (check)

INTEGER *4, INTENT (INOUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_multiplex_init (check)

INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_shutdown ()
END SUBROUTINE
SUBROUTINE PAPIF_create_eventset (EventSet, check)

INTEGER *4, INTENT (INOUT) :: EventSet
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_cleanup_eventset (EventSet, check)

INTEGER *4, INTENT (IN) :: EventSet
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_destroy_eventset (EventSet, check)

INTEGER *4, INTENT (INOUT) :: EventSet
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_set_multiplex (EventSet)

INTEGER *4, INTENT (IN) :: EventSet
END SUBROUTINE
SUBROUTINE PAPIF_add_event (EventSet, EventCode, check)

INTEGER *4, INTENT (IN) :: EventSet
INTEGER *4, INTENT (IN) :: EventCode
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_list_events (EventSet, Events, number, check)

INTEGER *4, INTENT (IN) :: EventSet
INTEGER *4, INTENT (INOUT) :: number
INTEGER *4, DIMENSION (number), INTENT (OUT) :: Events
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_start (EventSet, check)

INTEGER *4, INTENT (IN) :: EventSet
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_reset (EventSet, check)

INTEGER *4, INTENT (IN) :: EventSet
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_stop (EventSet, values, check)

INTEGER *4, INTENT (IN) :: EventSet
INTEGER *8, DIMENSION (*), INTENT (OUT) :: values
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_read (EventSet, values, check)

INTEGER *4, INTENT (IN) :: EventSet
INTEGER *8, DIMENSION (*), INTENT (OUT) :: values
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_accum (EventSet, values, check)

INTEGER *4, INTENT (IN) :: EventSet
INTEGER *8, DIMENSION (*), INTENT (OUT) :: values
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_get_clockrate (clockrate)

INTEGER *4, INTENT (OUT) :: clockrate
END SUBROUTINE
SUBROUTINE PAPIF_get_real_cyc (real_cyc)

INTEGER *8 :: real_cyc
END SUBROUTINE
SUBROUTINE PAPIF_get_virt_cyc (virt_cyc)

INTEGER *8 :: virt_cyc
END SUBROUTINE
SUBROUTINE PAPIF_event_code_to_name (EventCode, EventName, check)

INTEGER *4, INTENT (IN) :: EventCode
CHARACTER ( LEN=*), INTENT (OUT) :: EventName
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_event_name_to_code (EventName, EventCode, check)

CHARACTER ( LEN=*), INTENT (IN) :: EventName
INTEGER *4, INTENT (OUT) :: EventCode
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
SUBROUTINE PAPIF_perror (code, destination, check)

INTEGER *4, INTENT (IN) :: code
CHARACTER ( LEN=*), INTENT (OUT) :: destination
INTEGER *4, INTENT (OUT) :: check
END SUBROUTINE
