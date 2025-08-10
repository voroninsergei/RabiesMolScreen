def test_help():
 import subprocess,sys
 r=subprocess.run([sys.executable,'-m','rabiesmol.cli','--help'])
 assert r.returncode==0
