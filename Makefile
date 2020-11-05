Main_Name := main
MyGF2m := myGF2m

App_C_Flags := -g -O0 -Wall -Wextra -Wvla -Wno-unknown-pragmas -I.
App_Cpp_Flags := $(App_C_Flags) -std=c++14
App_Link_Flags := -lcrypto

all: $(Main_Name)

mygf2m: $(MyGF2m)

reed_solomon.o: reed_solomon.c reed_solomon.h 
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

$(Main_Name): reed_solomon.o
	@$(CXX) $^ -o $@ $(App_Link_Flags)
	@echo "LINK =>  $@"

myGF2m.o: myGF2m.c myGF2m.h
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

reed_solomon_myssl.o: reed_solomon_myssl.c reed_solomon_myssl.h myGF2m.o
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

$(MyGF2m): myGF2m.o reed_solomon_myssl.o 
	@$(CXX) $^ -o $@ $(App_Link_Flags)
	@echo "LINK =>  $@"

clean:
	@rm -rf $(Main_Name) *.o