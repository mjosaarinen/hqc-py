test:
	grep `python3 hqc.py | sha256sum | colrm 65` kat/*

clean:
	$(RM) -f *.pyc *.cprof */*.pyc *.rsp *.log
	$(RM) -rf __pycache__ */__pycache__

