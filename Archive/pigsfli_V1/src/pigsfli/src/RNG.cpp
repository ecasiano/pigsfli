#include "RNG.h"

//DEFINATIONS OF FUNCTIONS FROM MTfromPIMC CLASS

	 MTFromPIMC::MTFromPIMC(){ seed();};
	
	
   	 MTFromPIMC::MTFromPIMC(const uint32 _seed){seed (_seed);};

	 void MTFromPIMC::seed(){
	random.seed();
	};
	
     void MTFromPIMC::seed(const uint32 seed){
	random.seed(seed);
	};
	 uint32 MTFromPIMC::randInt(){return random.randInt();};
	 uint32 MTFromPIMC::randInt( const uint32 n ){return random.randInt( n );};
	 double MTFromPIMC::rand(){return random.rand();};
	
	
	 std::stringstream MTFromPIMC::save(){
	    std::stringstream stateStrStrm;
        stateStrStrm.str("");

        /* Save the state of the random number generator */
        uint32 randomState[random.SAVE];
        random.save(randomState);
        for (int i = 0; i < random.SAVE; i++)
            stateStrStrm << randomState[i] << " ";
        stateStrStrm << std::endl;
    
        return stateStrStrm ;
	};
	
	
	 void MTFromPIMC::load(std::istream& stateStrStrm){
	uint32 randomState[random.SAVE];
    for (int i = 0; i < random.SAVE; i++) 
        stateStrStrm >> randomState[i];
    random.load(randomState);
	};
	
     double MTFromPIMC::randExc(){
    return random.randExc(); 
    };
    
     double MTFromPIMC::randExc( const double n )
    { return random.randExc(n); };

     double MTFromPIMC::randDblExc()
    { return  random.randDblExc(); };

     double MTFromPIMC::randDblExc( const double n )
    { return random.randDblExc(n); };
    
     double MTFromPIMC::randNorm( const double mean, const double stddev = 1.0 ){return random.randNorm(mean , stddev);};
     
     
     
     
// DEFINITIONS OF FUNCTIONS FROM MTFromSTL

	 MTFromSTL::MTFromSTL(){
	    seed();
	    using parm_t = decltype(dis)::param_type;
	    dis.param(parm_t{0.0, 1.0});
	};
     MTFromSTL::MTFromSTL(const uint32 _seed){
        seed (_seed);
        using parm_t = decltype(dis)::param_type;
	    dis.param(parm_t{0.0, 1.0});
    };
    

	
	 void MTFromSTL::seed(){generator.seed();};
	 void MTFromSTL::seed( const uint32 seed ){generator.seed(seed);};

	 uint32 MTFromSTL::randInt(){
	    return generator();
	};
	 uint32 MTFromSTL::randInt( const uint32 n ){
	    using parm_t = decltype(disInt)::param_type;
	    return disInt( generator, parm_t{0, n});
	};
	
	
	
	
	 double MTFromSTL::rand(){return randExc();};
	 
	 
	 
	 std::stringstream MTFromSTL::save(){
	    std::stringstream stateStrStrm;
        stateStrStrm.str("");

        /* Save the state of the random number generator */
       
        stateStrStrm << generator;
        stateStrStrm << std::endl;
        stateStrStrm << disInt;
        stateStrStrm << std::endl;
        stateStrStrm << dis;
        stateStrStrm << std::endl;
        stateStrStrm << disNorm;
        stateStrStrm << std::endl;   
        return stateStrStrm ;
	};
	
	
	 void MTFromSTL::load(std::istream& stateStrStrm){
	    stateStrStrm >> generator;
        stateStrStrm >> disInt;
        stateStrStrm >> dis;
        stateStrStrm >> disNorm;   

	};
	 
	 double MTFromSTL::randExc(){
    return dis(generator); 
    };


    double MTFromSTL::randExc( const double n ){
	    using parm_t = decltype(dis)::param_type;
	    return dis( generator, parm_t{0, n});
	};

    double MTFromSTL::randNorm( const double mean, const double stddev = 1.0){
	    using parm_t = decltype(disNorm)::param_type;
	    return disNorm( generator, parm_t{mean, stddev});
	};

     double MTFromSTL::randDblExc()
    { return  randExc(); };

     double MTFromSTL::randDblExc( const double n )
    { return randExc(n); };
     
     
// DEFINITIONS OF FUNCTIONS FROM MTFromPCG

	 MTFromPCG::MTFromPCG(){
	    seed();
	    using parm_t = decltype(dis)::param_type;
	    dis.param(parm_t{0.0, 1.0});
	};
     MTFromPCG::MTFromPCG(const uint32 _seed){
        seed (_seed);
        using parm_t = decltype(dis)::param_type;
	    dis.param(parm_t{0.0, 1.0});
    };
    

	
	 void MTFromPCG::seed(){generator.seed();};
	 void MTFromPCG::seed( const uint32 seed ){generator.seed(seed);};

	 uint32 MTFromPCG::randInt(){
	    return generator();
	};
	 uint32 MTFromPCG::randInt( const uint32 n ){
	    using parm_t = decltype(disInt)::param_type;
	    return disInt( generator, parm_t{0, n});
	};
	
	
	
	
	 double MTFromPCG::rand(){return randExc();};
	 
	 
	 
	 std::stringstream MTFromPCG::save(){
	    std::stringstream stateStrStrm;
        stateStrStrm.str("");

        /* Save the state of the random number generator */
       
        stateStrStrm << generator;
        stateStrStrm << std::endl;
        stateStrStrm << disInt;
        stateStrStrm << std::endl;
        stateStrStrm << dis;
        stateStrStrm << std::endl;
        stateStrStrm << disNorm;
        stateStrStrm << std::endl;   
        return stateStrStrm ;
	};
	
	
	 void MTFromPCG::load(std::istream& stateStrStrm){
	    stateStrStrm >> generator;
        stateStrStrm >> disInt;
        stateStrStrm >> dis;
        stateStrStrm >> disNorm;   

	};
	 
	 double MTFromPCG::randExc(){
    return dis(generator); 
    };


    double MTFromPCG::randExc( const double n ){
	    using parm_t = decltype(dis)::param_type;
	    return dis( generator, parm_t{0, n});
	};

    double MTFromPCG::randNorm( const double mean, const double stddev = 1.0){
	    using parm_t = decltype(disNorm)::param_type;
	    return disNorm( generator, parm_t{mean, stddev});
	};

     double MTFromPCG::randDblExc()
    { return  randExc(); };

     double MTFromPCG::randDblExc( const double n )
    { return randExc(n); };
     
     
// DEFINITIONS OF FUNCTIONS FROM MTFromBOOST

	 MTFromBOOST::MTFromBOOST(){
	    seed();
	    using parm_t = decltype(dis)::param_type;
	    dis.param(parm_t{0.0, 1.0});
	};
     MTFromBOOST::MTFromBOOST(const uint32 _seed){
        seed (_seed);
        using parm_t = decltype(dis)::param_type;
	    dis.param(parm_t{0.0, 1.0});
    };
    

	
	 void MTFromBOOST::seed(){generator.seed();};
	 void MTFromBOOST::seed( const uint32 seed ){generator.seed(seed);};

	 uint32 MTFromBOOST::randInt(){
	    return generator();
	};
	 uint32 MTFromBOOST::randInt( const uint32 n ){
	    using parm_t = decltype(disInt)::param_type;
	    return disInt( generator, parm_t{0, n});
	};
	
	
	 double MTFromBOOST::rand(){return randExc();};
	 
	 std::stringstream MTFromBOOST::save(){
	    std::stringstream stateStrStrm;
        stateStrStrm.str("");

        /* Save the state of the random number generator */
       
        stateStrStrm << generator;
        stateStrStrm << std::endl;
        stateStrStrm << disInt;
        stateStrStrm << std::endl;
        stateStrStrm << dis;
        stateStrStrm << std::endl;
        stateStrStrm << disNorm;
        stateStrStrm << std::endl;    
        return stateStrStrm ;
	};
	
	
	 void MTFromBOOST::load(std::istream& stateStrStrm){
	    stateStrStrm >> generator;
        stateStrStrm >> disInt;
        stateStrStrm >> dis;    
        stateStrStrm >> disNorm;

	};
	 
	 double MTFromBOOST::randExc(){
    return dis(generator); 
    };


    double MTFromBOOST::randExc( const double n ){
	    using parm_t = decltype(dis)::param_type;
	    return dis( generator, parm_t{0, n});
	};


     double MTFromBOOST::randDblExc()
    { return  randExc(); };

     double MTFromBOOST::randDblExc( const double n )
    { return randExc(n); };
    	
     double MTFromBOOST::randNorm( const double mean, const double stddev = 1.0){
	    using parm_t = decltype(disNorm)::param_type;
	    return disNorm( generator, parm_t{mean, stddev});
	};
