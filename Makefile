SOURCE_DIR = src
LINK_DIR = obj

.PHONY: all clean

all : rf_dist rfa_dist l1_dist l2_dist quartet_dist

clean :
	rm -rf $(LINK_DIR) rf_dist rfa_dist l1_dist l2_dist quartet_dist

$(SOURCE_DIR) :
	if ! [ -d $(SOURCE_DIR) ]; then echo "Source dir doesn't exist"; exit 1; fi

$(LINK_DIR) :
	if ! [ -d $(LINK_DIR) ]; then mkdir $(LINK_DIR); fi

treedist.o : $(SOURCE_DIR) $(LINK_DIR) $(SOURCE_DIR)/treedist.c
	gcc -O2 -c $(SOURCE_DIR)/treedist.c -o $(LINK_DIR)/treedist.o

rf_dist.o : $(SOURCE_DIR) $(LINK_DIR) $(SOURCE_DIR)/rf_dist.c
	gcc -O2 -c $(SOURCE_DIR)/rf_dist.c -o $(LINK_DIR)/rf_dist.o

rfa_dist.o : $(SOURCE_DIR) $(LINK_DIR) $(SOURCE_DIR)/rfa_dist.c
	gcc -O2 -c $(SOURCE_DIR)/rfa_dist.c -o $(LINK_DIR)/rfa_dist.o

l1_dist.o : $(SOURCE_DIR) $(LINK_DIR) $(SOURCE_DIR)/l1_dist.c
	gcc -O2 -c $(SOURCE_DIR)/l1_dist.c -o $(LINK_DIR)/l1_dist.o

l2_dist.o : $(SOURCE_DIR) $(LINK_DIR) $(SOURCE_DIR)/l2_dist.c
	gcc -O2 -c $(SOURCE_DIR)/l2_dist.c -o $(LINK_DIR)/l2_dist.o

quartet_dist.o : $(SOURCE_DIR) $(LINK_DIR) $(SOURCE_DIR)/quartet_dist.c
	gcc -O2 -c $(SOURCE_DIR)/quartet_dist.c -o $(LINK_DIR)/quartet_dist.o

rf_dist : rf_dist.o treedist.o
	gcc $(LINK_DIR)/rf_dist.o $(LINK_DIR)/treedist.o -lm -o rf_dist
 
rfa_dist : rfa_dist.o treedist.o
	gcc $(LINK_DIR)/rfa_dist.o $(LINK_DIR)/treedist.o -lm -o rfa_dist

l1_dist : l1_dist.o treedist.o
	gcc $(LINK_DIR)/l1_dist.o $(LINK_DIR)/treedist.o -lm -o l1_dist

l2_dist : l2_dist.o treedist.o
	gcc $(LINK_DIR)/l2_dist.o $(LINK_DIR)/treedist.o -lm -o l2_dist

quartet_dist : quartet_dist.o treedist.o
	gcc $(LINK_DIR)/quartet_dist.o $(LINK_DIR)/treedist.o -lm -o quartet_dist
