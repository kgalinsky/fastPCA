MAKE_PATH=$(dir $(lastword $(MAKEFILE_LIST)))
include $(MAKE_PATH)/vars.mk

PROJ_PATH=$(MAKE_PATH)/..
export PATH:=$(PROJ_PATH)/src:$(PROJ_PATH)/scripts:$(PROJ_PATH)/R:$(PATH)

.SECONDARY:

%.pca.evec %.pca.eval: $$(b1).geno
	pca -o $*.pca -k $(s1) $<

%.fast.evec: $$(b4).geno
	GSL_RNG_SEED=$(s1) fastpca -k $(s4) -l $(s3) -i $(s2) -o $*.fast $<

ladder/%.evec: $$(b4).geno
	GSL_RNG_SEED=$(s2) fastpca-iterative -k $(s4) -l $(s3) -i $(s1) -o ladder/$(b1) $<

# One population; sim1.$(SNP).$(IND).geno
sim1.%.geno:
	rand_pop.pl $(t1) $(t2) > $@ 2> /dev/null

# Two populations of size N; sim2.$(SNP).$(FST).$(N).geno
sim2.%.geno:
	rand_pop.pl $(t1) [.$(t2),$(t3),$(t3)] > $@ 2> /dev/null

# Three populations of size N from one ancestor
sim3.%.geno:
	rand_pop.pl $(t1) [.$(t2),$(t3),$(t3),$(t3)] > $@ 2> /dev/null

# Two splits; sim3a.$(SNP).$(FST1).$(FST2).$(N).geno
sim3a.%.geno:
	rand_pop.pl $(t1) [.$(t2),[.$(t3),$(t4),$(t4)],$(t4)] > $@ 2> /dev/null

# Three splits
sim4a.%.geno:
	rand_pop.pl $(t1) [.$(t2),[.$(t3),[.$(t4),$(t5),$(t5)],$(t5)],$(t5)] > $@ 2> /dev/null

%.r2.trace: $$(b3).10.pca.evec ladder/%.10.evec
	i=0; while [ $$i -le 10 ] ; do \
		j=$$(printf %02d $$i) ; \
		echo -ne $$i"\t"; \
		r2.R $(a1) ladder/$*.$$j.evec | paste -s ; \
		(( i = i + 1 )) ; \
	done > $@

%.FoMKG.trace: $$(b3).10.pca.evec ladder/%.10.evec
	i=0; while [ $$i -le 10 ] ; do \
		j=$$(printf %02d $$i) ; \
		echo -ne $$i"\t"; \
		FoMKG.R $(a1) ladder/$*.$$j.evec | paste -s ; \
		(( i = i + 1 )) ; \
	done > $@

%.multitrace: $$(b1).000000.$$(s1).trace \
         $$(b1).100000.$$(s1).trace \
         $$(b1).200000.$$(s1).trace \
         $$(b1).300000.$$(s1).trace \
         $$(b1).400000.$$(s1).trace \
         $$(b1).500000.$$(s1).trace \
         $$(b1).600000.$$(s1).trace \
         $$(b1).700000.$$(s1).trace \
         $$(b1).800000.$$(s1).trace \
         $$(b1).900000.$$(s1).trace
	for i in $^; do echo $$i; cat $$i; echo; done > $@
