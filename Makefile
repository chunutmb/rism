default:
	$(MAKE) -C prog

clean:
	rm -rf $(bins) $(bins_d) *~ a.out *.dat* gmon.out \
	  .*.un~ .*~ */.*.un~ */*~ r[0-9]*hs
	$(MAKE) -C web clean
	$(MAKE) -C prog clean
	rstrip.py -Rlv

excludes = --exclude=".*" --exclude="*~" --exclude="bak"

Bossman: clean
	rsync -avzL $(excludes) * /Bossman/cz1/rism/

Bossman2: clean
	rsync -vzL $(excludes) * cz1@129.109.88.204:/Bossman/cz1/rism/


