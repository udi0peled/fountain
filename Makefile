BIGNUM_Name := bignum
MyGF2m_Name := myGF2m
Galois_Name := galois

App_C_Flags := -O2 -Wall -Wextra -Wvla -Wno-unknown-pragmas -I.
App_Cpp_Flags := $(App_C_Flags) -std=c++14  -march=skylake
App_Link_Flags := -lcrypto

ALL = bn mygf2m galois

bn	  : $(BIGNUM_Name)
mygf2m: $(MyGF2m_Name)
galois: $(Galois)

clean:
	@rm -rf $(ALL) *.o *.out

reed_solomon.o: reed_solomon.c reed_solomon.h 
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

$(BIGNUM_Name): reed_solomon.o
	@$(CXX) $^ -o $@.out $(App_Link_Flags)
	@echo "LINK =>  $@"

myGF2m.o: myGF2m.c myGF2m.h
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

reed_solomon_myGf2m.o: reed_solomon_myGF2m.c reed_solomon_myGF2m.h myGF2m.o
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

$(MyGF2m_Name): myGF2m.o reed_solomon_myGF2m.o 
	@$(CXX) $^ -o $@.out
	@echo "LINK =>  $@"

galois.o: galois.c galois.h
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

reed_solomon_galois.o: reed_solomon_galois.c reed_solomon_galois.h galois.o
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

$(Galois): galois.o reed_solomon_galois.o 
	@$(CXX) $^ -o $@.out
	@echo "LINK =>  $@"