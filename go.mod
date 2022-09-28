module git.sr.ht/~vejnar/GeneAbacus

go 1.19

//replace git.sr.ht/~vejnar/GeneAbacus => ./GeneAbacus

require (
	github.com/biogo/hts v1.4.4
	github.com/biogo/store v0.0.0-20201120204734-aad293a2328f
	github.com/klauspost/compress v1.15.11
	github.com/pierrec/lz4 v2.6.1+incompatible
	golang.org/x/sync v0.0.0-20220923202941-7f9b1623fab7
	gopkg.in/fatih/set.v0 v0.2.1
)

require github.com/frankban/quicktest v1.14.3 // indirect
