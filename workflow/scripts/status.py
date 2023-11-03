#!/usr/bin/env python
import subprocess
import sys

jobid = sys.argv[1]

output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

running_status=[ "PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED", "STOPPED", "SIGNALING", "RESIZING", "REQUEUED", "REQUEUE_HOLD", "REQUEUE_FED", "RESV_DEL_HOLD", "STAGE_OUT" ]

failure_states=[ "FAILED", "PREEMPTED", "TIMEOUT", "SPECIAL_EXIT", "REVOKED", "OUT_OF_MEMORY", "NODE_FAIL", "DEADLINE", "CANCELLED", "BOOT_FAIL" ]

if "COMPLETED" in output:
    print("success")
elif any(r in output for r in running_status):
    print("running")
elif any(r in output for r in failure_states):
    print("failed")
else:
    print("running")

