# build package documentation
all: doc vignette
.PHONY: doc test vignette

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

doc:
	R -e 'devtools::document()'


test:
	R -e 'devtools::test()'

vignette: docs/xoi.html

docs/xoi.html: docs/xoi.Rmd
	cd $(<D);R $(R_OPTS) -e "rmarkdown::render('$(<F)')"
