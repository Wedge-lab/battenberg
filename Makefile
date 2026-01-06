.PHONY: style lint test deps check install docs pak

docs:
	Rscript -e "roxygen2::roxygenise(clean = TRUE, roclets = c('rd', 'namespace'))"

# Run the auto-formatter (styler)
style:
	Rscript -e "styler::style_pkg(transformers = styler::tidyverse_style(strict = TRUE))"

# Run the linter
lint:
	Rscript -e "lintr::lint_package()"

# Run the full test suite
test:
	Rscript -e "devtools::test()"

pak:
	@echo "Installing pak and core dependencies..."
	Rscript -e "install.packages('pak', repos='https://cloud.r-project.org')"

deps:
	@echo "Installing all dependencies listed in DESCRIPTION..."
	export GITHUB_PAT=""; \
	Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); \
		            pak::repo_add(Bioc = '3.18'); \
		            pak::local_install_deps(upgrade = TRUE, dependencies = TRUE)"
check:
	Rscript -e "devtools::check(error_on = 'warning')"

install:
	@echo "Installing Battenberg..."
	Rscript -e "remotes::install_local('.', upgrade='never')"
