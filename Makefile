BIGNUM_Name := bignum_rs
MyGF2m_Name := myGF2m_rs
Galois_Name := galois_rs

App_C_Flags := -g -O0 -Wall -Wextra -Wvla -Wno-unknown-pragmas -I.
App_Cpp_Flags := $(App_C_Flags) -std=c++14 
App_Link_Flags := -lcrypto

ALL := bn mygf2m galois

bn	 : $(BIGNUM_Name)
mygf : $(MyGF2m_Name)
gal  : $(Galois_Name)

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

reed_solomon_galois.o: reed_solomon_galois.c reed_solomon_galois.h galois_16bit_log_table.h  galois_16bit_inv_log_table.h 
	@$(CC) $(App_C_Flags) -c $< -o $@
	@echo "CC   <=  $<"

$(Galois_Name): reed_solomon_galois.o 
	@$(CXX) $^ -o $@.out
	@echo "LINK =>  $@"
