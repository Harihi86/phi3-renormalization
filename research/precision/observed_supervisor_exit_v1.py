"""Read the real supervisor exit after a child-exit/status-write race.

This observer sends no signals and never infers success from child absence.
"""
import time


def observe(path,expected_entry,identity,read,*,grace=3,clock=time.monotonic,sleep=time.sleep):
    deadline=clock()+grace
    while True:
        state=read(path)
        if state['state']=='terminal':
            assert state['exit_code']==0 and state['sources_unchanged'],'Observed unsuccessful supervisor exit'
            return state
        assert state['state']=='running'
        child=state['child_pid'];assert type(child) is int and child>0
        if identity(child).endswith(str(expected_entry)):return state
        assert clock()<deadline,'No successful supervisor receipt after child identity disappeared'
        sleep(.05)
