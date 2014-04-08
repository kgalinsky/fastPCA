# munge target to create variables

ifndef _VARS
_VARS 	= 1

# tokens
.SECONDEXPANSION:
tokens	= $(subst ., ,$*)
t1	= $(word 1,$(tokens))
t2	= $(word 2,$(tokens))
t3	= $(word 3,$(tokens))
t4	= $(word 4,$(tokens))
t5	= $(word 5,$(tokens))
t6	= $(word 6,$(tokens))
t7	= $(word 7,$(tokens))
t8	= $(word 8,$(tokens))
t9	= $(word 9,$(tokens))

# basename
.SECONDEXPANSION:
b1      = $(basename $*)
b2      = $(basename $(b1))
b3      = $(basename $(b2))
b4      = $(basename $(b3))
b5      = $(basename $(b4))
b6      = $(basename $(b5))
b7      = $(basename $(b6))
b8      = $(basename $(b7))
b9      = $(basename $(b8))

# suffix
.SECONDEXPANSION:
s1      = $(strip $(subst ., ,$(suffix $*)))
s2      = $(strip $(subst ., ,$(suffix $(b1))))
s3      = $(strip $(subst ., ,$(suffix $(b2))))
s4      = $(strip $(subst ., ,$(suffix $(b3))))
s5      = $(strip $(subst ., ,$(suffix $(b4))))
s6      = $(strip $(subst ., ,$(suffix $(b5))))
s7      = $(strip $(subst ., ,$(suffix $(b6))))
s8      = $(strip $(subst ., ,$(suffix $(b7))))
s9      = $(strip $(subst ., ,$(suffix $(b8))))

# arguments
.SECONDEXPANSION:
a1	= $(word 1,$^)
a2	= $(word 2,$^)
a3	= $(word 3,$^)
a4	= $(word 4,$^)
a5	= $(word 5,$^)
a6	= $(word 6,$^)
a7	= $(word 7,$^)
a8	= $(word 8,$^)
a9	= $(word 9,$^)

endif
