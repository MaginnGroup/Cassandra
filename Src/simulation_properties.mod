	  þ  G   k820309              14.0         7íT                                                                                                           
       simulation_properties.f90 SIMULATION_PROPERTIES                                                    
                                                          
                         @                                '                    #MOLECULE_TYPE    #RX_NUM    #LIVE    #INSIDE    #WHICH_BOX    #XCOM 	   #YCOM 
   #ZCOM    #EULER1    #EULER2    #EULER3    #XCOM_OLD    #YCOM_OLD    #ZCOM_OLD    #EULER1_OLD    #EULER2_OLD    #EULER3_OLD    #CFC_LAMBDA    #MAX_DCOM    #MAX_DCOM_OLD                                                                                                                                                                                                                                                                                                                                                                      	               
                                              
                
                                                   (          
                                                   0       	   
                                                   8       
   
                                                   @          
                                                   H          
                                                   P          
                                                   X          
                                                   `          
                                                   h          
                                                   p          
                                                   x          
                                                             
                                                             
                     @                                '                    #VDW_POTENTIAL_TYPE    #VDW_PARAM    #ELEMENT    #ATOM_NAME    #MASS    #CHARGE    #ATOM_TYPE_NUMBER    #RING_ATOM                                                                                                                           
                        
  p          p 
           p 
                                                                                h                                                                     j                                                              p          
                                                   x          
                                                                                                                                                                          !                                   &                                                                                     "                                   &                   &                                                                                     #                                    &                   &                                           #MOLECULE_CLASS                                                 $                                                       0         @                                %                                   &                   &                                                                                       &                                                      '                                   &                                                                                     (                                    &                   &                                           #NONBOND_CLASS             @@                               )                                   &                                                                                       *                     @@                               +                                   &                                                    @@                               ,                                   &                                           #         @                                   -                    #THIS_BOX .   #THIS_SPECIES /   #NMOLECULES_SPECIES 0                                              .                                                       /                      D                                 0            #         @                                   1                    #THIS_BOX 2   #THIS_SPECIES 3   #IM 4   #ALIVE 5                                              2                                                       3                                                       4                      D                                 5            #         @                                   6                    #THIS_BOX 7   #IS 8   #IM 9   #POSITION :                                              7                                                       8                                                       9                      D                                 :            #         @                                   ;                    #THIS_BOX <   #THIS_SPECIES =   #IM >   #ALIVE ?                                              <                                                       =                                                       >                      D                                 ?            #         @                                   @                    #IS A   #ALIVE B   #POSITION C                                              A                                                       B                      D                                 C            #         @                                   D                   #COMPUTE_BEADS%ALLOCATED E   #THIS_BOX F                                              E     ALLOCATED           
                                  F                  8      fn#fn !   Ø   @   J   TYPE_DEFINITIONS      @   J   RUN_VARIABLES 0   X  `      MOLECULE_CLASS+TYPE_DEFINITIONS >   ¸  H   a   MOLECULE_CLASS%MOLECULE_TYPE+TYPE_DEFINITIONS 7      H   a   MOLECULE_CLASS%RX_NUM+TYPE_DEFINITIONS 5   H  H   a   MOLECULE_CLASS%LIVE+TYPE_DEFINITIONS 7     H   a   MOLECULE_CLASS%INSIDE+TYPE_DEFINITIONS :   Ø  H   a   MOLECULE_CLASS%WHICH_BOX+TYPE_DEFINITIONS 5      H   a   MOLECULE_CLASS%XCOM+TYPE_DEFINITIONS 5   h  H   a   MOLECULE_CLASS%YCOM+TYPE_DEFINITIONS 5   °  H   a   MOLECULE_CLASS%ZCOM+TYPE_DEFINITIONS 7   ø  H   a   MOLECULE_CLASS%EULER1+TYPE_DEFINITIONS 7   @  H   a   MOLECULE_CLASS%EULER2+TYPE_DEFINITIONS 7     H   a   MOLECULE_CLASS%EULER3+TYPE_DEFINITIONS 9   Ð  H   a   MOLECULE_CLASS%XCOM_OLD+TYPE_DEFINITIONS 9     H   a   MOLECULE_CLASS%YCOM_OLD+TYPE_DEFINITIONS 9   `  H   a   MOLECULE_CLASS%ZCOM_OLD+TYPE_DEFINITIONS ;   ¨  H   a   MOLECULE_CLASS%EULER1_OLD+TYPE_DEFINITIONS ;   ð  H   a   MOLECULE_CLASS%EULER2_OLD+TYPE_DEFINITIONS ;   8  H   a   MOLECULE_CLASS%EULER3_OLD+TYPE_DEFINITIONS ;     H   a   MOLECULE_CLASS%CFC_LAMBDA+TYPE_DEFINITIONS 9   È  H   a   MOLECULE_CLASS%MAX_DCOM+TYPE_DEFINITIONS =     H   a   MOLECULE_CLASS%MAX_DCOM_OLD+TYPE_DEFINITIONS /   X  Î       NONBOND_CLASS+TYPE_DEFINITIONS B   &	  P   a   NONBOND_CLASS%VDW_POTENTIAL_TYPE+TYPE_DEFINITIONS 9   v	     a   NONBOND_CLASS%VDW_PARAM+TYPE_DEFINITIONS 7   
  P   a   NONBOND_CLASS%ELEMENT+TYPE_DEFINITIONS 9   b
  P   a   NONBOND_CLASS%ATOM_NAME+TYPE_DEFINITIONS 4   ²
  H   a   NONBOND_CLASS%MASS+TYPE_DEFINITIONS 6   ú
  H   a   NONBOND_CLASS%CHARGE+TYPE_DEFINITIONS @   B  H   a   NONBOND_CLASS%ATOM_TYPE_NUMBER+TYPE_DEFINITIONS 9     H   a   NONBOND_CLASS%RING_ATOM+TYPE_DEFINITIONS )   Ò         NMOLECULES+RUN_VARIABLES %   ^  ¤       LOCATE+RUN_VARIABLES ,     ¸       MOLECULE_LIST+RUN_VARIABLES )   º  q       INT_NORMAL+RUN_VARIABLES )   +  ¤       NINT_BEADS+RUN_VARIABLES '   Ï  @       NSPECIES+RUN_VARIABLES %            NATOMS+RUN_VARIABLES +     ·       NONBOND_LIST+RUN_VARIABLES (   R         NBEADS_IN+RUN_VARIABLES ,   Þ  @       NBR_ATOMTYPES+RUN_VARIABLES )            NBEADS_OUT+RUN_VARIABLES ,   ª         NBEADSFRAC_IN+RUN_VARIABLES '   6         GET_NMOLECULES_SPECIES 0   ¶  @   a   GET_NMOLECULES_SPECIES%THIS_BOX 4   ö  @   a   GET_NMOLECULES_SPECIES%THIS_SPECIES :   6  @   a   GET_NMOLECULES_SPECIES%NMOLECULES_SPECIES #   v  {       GET_INDEX_MOLECULE ,   ñ  @   a   GET_INDEX_MOLECULE%THIS_BOX 0   1  @   a   GET_INDEX_MOLECULE%THIS_SPECIES &   q  @   a   GET_INDEX_MOLECULE%IM )   ±  @   a   GET_INDEX_MOLECULE%ALIVE &   ñ  t       GET_POSITION_MOLECULE /   e  @   a   GET_POSITION_MOLECULE%THIS_BOX )   ¥  @   a   GET_POSITION_MOLECULE%IS )   å  @   a   GET_POSITION_MOLECULE%IM /   %  @   a   GET_POSITION_MOLECULE%POSITION +   e  {       GET_INDEX_INTEGER_MOLECULE 4   à  @   a   GET_INDEX_INTEGER_MOLECULE%THIS_BOX 8      @   a   GET_INDEX_INTEGER_MOLECULE%THIS_SPECIES .   `  @   a   GET_INDEX_INTEGER_MOLECULE%IM 1      @   a   GET_INDEX_INTEGER_MOLECULE%ALIVE $   à  i       GET_POSITION_LOCATE '   I  @   a   GET_POSITION_LOCATE%IS *     @   a   GET_POSITION_LOCATE%ALIVE -   É  @   a   GET_POSITION_LOCATE%POSITION    	  s       COMPUTE_BEADS (   |  B      COMPUTE_BEADS%ALLOCATED '   ¾  @   a   COMPUTE_BEADS%THIS_BOX 