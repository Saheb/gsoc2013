all: all-how-tos all-regression-tests

clean: clean-how-tos clean-regression-tests

all-%:
	$(MAKE) -C $(patsubst all-%,%,$@) all

clean-%:
	$(MAKE) -C $(patsubst clean-%,%,$@) clean
