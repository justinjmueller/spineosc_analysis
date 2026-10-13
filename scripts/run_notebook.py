#!/usr/bin/env python3
"""Executes a notebook's code cells in a Jupyter kernel and writes the outputs back (no nbclient/nbconvert needed).
Usage (mutools env):  python scripts/run_notebook.py plotting/systematics.ipynb [kernel-name] [--timeout S]"""
import queue, sys, time
import nbformat
from jupyter_client.manager import KernelManager

args = sys.argv[1:]; timeout = 3600
if "--timeout" in args: i = args.index("--timeout"); timeout = int(args[i + 1]); del args[i:i + 2]
path = args[0]; kernel = args[1] if len(args) > 1 else "mutools-dev"
nb = nbformat.read(path, as_version=4)
km = KernelManager(kernel_name=kernel); km.start_kernel(cwd=str(__import__("pathlib").Path(path).resolve().parent))
kc = km.client(); kc.start_channels(); kc.wait_for_ready(timeout=120)
failed = 0
for i, cell in enumerate(nb.cells):
    if cell.cell_type != "code" or not cell.source.strip():
        continue
    cell.outputs = []; msg_id = kc.execute(cell.source); t0 = time.time(); done = False
    while not done:
        try: msg = kc.get_iopub_msg(timeout=timeout)
        except queue.Empty: print(f"cell {i}: timed out"); failed += 1; break
        if msg["parent_header"].get("msg_id") != msg_id: continue
        t, c = msg["msg_type"], msg["content"]
        if t == "stream": cell.outputs.append(nbformat.v4.new_output("stream", name=c["name"], text=c["text"]))
        elif t in ("display_data", "execute_result"):
            cell.outputs.append(nbformat.v4.new_output(t, data=c["data"], metadata=c.get("metadata", {}), **({"execution_count": c.get("execution_count")} if t == "execute_result" else {})))
        elif t == "error": cell.outputs.append(nbformat.v4.new_output("error", ename=c["ename"], evalue=c["evalue"], traceback=c["traceback"])); failed += 1
        elif t == "status" and c["execution_state"] == "idle": done = True
    print(f"cell {i}: {time.time()-t0:.0f}s, {len(cell.outputs)} outputs")
kc.stop_channels(); km.shutdown_kernel(now=True)
nbformat.write(nb, path); print("wrote", path, "| failed cells:", failed)
sys.exit(1 if failed else 0)
