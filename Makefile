SRC_DIR=src
MAIN=convolution

compile:
	mpic++ $(SRC_DIR)/$(MAIN).cpp -o $(MAIN).out
