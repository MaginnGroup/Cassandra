	  èW  Þ   k820309              14.0        ¼óT                                                                                                           
       fragment_growth.f90 FRAGMENT_GROWTH                                                   
                                                          
                                                          
                                                          
                         @                                'P              
      #RXP    #RYP    #RZP    #RXP_NLS 	   #RYP_NLS 
   #RZP_NLS    #RXP_OLD    #RYP_OLD    #RZP_OLD    #EXIST                                                               
                                                             
                                                             
                                              	               
                                              
                
                                                   (          
                                                   0          
                                                   8          
                                                   @       	   
                                                    H       
                        @                                '                    #MOLECULE_TYPE    #RX_NUM    #LIVE    #INSIDE    #WHICH_BOX    #XCOM    #YCOM    #ZCOM    #EULER1    #EULER2    #EULER3    #XCOM_OLD    #YCOM_OLD    #ZCOM_OLD    #EULER1_OLD    #EULER2_OLD     #EULER3_OLD !   #CFC_LAMBDA "   #MAX_DCOM #   #MAX_DCOM_OLD $                                                                                                                                                                                                                                                                                                                                                                                    
                                                              
                                                   (          
                                                   0       	   
                                                   8       
   
                                                   @          
                                                   H          
                                                   P          
                                                   X          
                                                   `          
                                                    h          
                                              !     p          
                                              "     x          
                                              #               
                                              $               
                     @                           %     '¨                   #BOX_SHAPE &   #INT_BOX_SHAPE '   #LENGTH (   #LENGTH_INV )   #MAX_DELTA *   #HLENGTH +   #BASIS_LENGTH ,   #COS_ANGLE -   #ANGLE .   #FACE_DISTANCE /   #VOLUME 0   #DV_MAX 1                                              &                                                                       '                                                             (     	                        
  p          p          p            p          p                                                                     )     	       `                 
  p          p          p            p          p                                                                     *     	       ¨                 
  p          p          p            p          p                                                                     +     	       ð                 
  p          p          p            p          p                                                                     ,            8                
  p          p            p                                                                     -            P                
  p          p            p                                                                     .            h             	   
  p          p            p                                                                     /                         
   
  p          p            p                                                                     0              
                                              1               
                     @               @           2     '                   #NATOMS 3   #NCONNECT 4   #NANCHORS 5   #TYPE 6   #NCONFIG 7   #ANCHOR 8   #ATOMS 9   #FRAG_CONNECT :   #RING ;   #RCUT_VDWSQ <   #RCUT_COULSQ =   #ALPHA_EWALD >                                               3                                                               4                                                              5                                                              6                                                              7                                                            8                                         &                                                                                     9            `                             &                                                                                     :            ¨                             &                                                                                       ;     ð       	                                                 <     ø       
   
                                              =               
                                              >              
                     @                           ?     '                    #FRAGMENT1 @   #FRAGMENT2 A                                               @                                                               A                                    @                           B     '                    #RXP C   #RYP D   #RZP E                                              C                
                                              D               
                                              E               
                                             F                                   &                                                                                        G                                                                     @                               H                        @                               I                        @                               J                      @                               K                                   &                                                      @                                L                     @                                M                                    &                   &                                           #MOLECULE_CLASS             @                                N            P                        &                   &                   &                                           #ATOM_CLASS                                                O            %         @                               P                    
       #RRANF%IAND Q   #RRANF%ISHFT R   #RRANF%IEOR S                                               Q     IAND                                             R     ISHFT                                             S     IEOR                                             T                                                      U                                   &                   &                                           #FRAG_CLASS 2                                             V                                    &                   &                   &                                           #FRAG_COORD_CLASS B                                               W                                                        X                                                      Y            ¨                       &                                           #BOX_CLASS %                                                Z                                                       0#         @                                  [                   #COMPUTE_MOLECULE_NONBOND_INTER_ENERGY%SUM \   #IM ]   #IS ^   #E_INTER_VDW _   #E_INTER_QQ `   #OVERLAP a                                               \     SUM           
                                  ]                     
                                  ^                                                     _     
                                                 `     
                                                  a                                                     b                   
                &                                                                                       c     
                
                      A@        35.0#         @                                  d                    #THIS_FRAG e   #THIS_IM f   #IS g   #THIS_BOX h   #NRG_RING_FRAG i             
                                  e                     
                                  f                     
                                  g                                                      h                                                      i     
                                                j                   
                &                   &                                           #         @                                  k                   #COMPUTE_MOLECULE_DIHEDRAL_ENERGY%COS l   #COMPUTE_MOLECULE_DIHEDRAL_ENERGY%DCOS m   #MOLECULE n   #SPECIES o   #ENERGY_DIHED p                                               l     COS                                             m     DCOS                                            n                                                       o                                                      p     
                                                 q                                   &                                                                                     r                                    &                   &                                           #FRAGMENT_BOND_CLASS ?   +           @@                              s     
       P             p          p 
           p 
                                                                             t     
                
                 àDTû!@        6.2831853072#         @                                  u                   #COMPUTE_ATOM_NONBOND_ENERGY%SQRT v   #COMPUTE_ATOM_NONBOND_ENERGY%EXP w   #THIS_ATOM x   #THIS_MOLECULE y   #THIS_SPECIES z   #E_INTRA_VDW {   #E_INTER_VDW |   #E_INTRA_QQ }   #E_INTER_QQ ~   #OVERLAP                                                v     SQRT                                             w     EXP           
                                  x                     
                                  y                     
                                  z                                                     {     
                                                 |     
                                                 }     
                                                 ~     
                                                              #         @                                                   	   #BUILD_MOLECULE%DEXP    #BUILD_MOLECULE%INT    #BUILD_MOLECULE%ALLOCATED    #BUILD_MOLECULE%MAXVAL    #BUILD_MOLECULE%MAX    #THIS_IM    #IS    #THIS_BOX    #FRAG_ORDER    #THIS_LAMBDA    #WHICH_ANCHOR    #ATTEMPT_PROB    #NRG_RING_FRAG_TOTAL    #CBMC_OVERLAP                                                    DEXP                                                 INT                                                 ALLOCATED                                                 MAXVAL                                                 MAX           D @                                                     D @                                                     D @                                                    D @                                                        p          & p        5 2 r F     5  p        r          5 2 r F     5  p        r                                  
  @                                   
                                                                       D @                                   
                 D @                                   
                 
D @                                           #         @                                                    #FRAGMENT_ORDER%INT    #FRAGMENT_ORDER%ALLOCATED    #THIS_FRAG    #IS    #FRAG_TOTAL    #FRAG_ORDER    #LIVE    #DEADEND    #CENTRAL                                                    INT                                                 ALLOCATED           D @                                                     D @                                                     D @                                                    D @                                                        p          & p        5 2 r F     5  p        r          5 2 r F     5  p        r                                 D @                                                        p          & p        5 2 r F     5  p        r          5 2 r F     5  p        r                                 D @                                                        p          & p        5 2 r F     5  p        r          5 2 r F     5  p        r                                 D @                                                        p          & p        5 2 r F     5  p        r          5 2 r F     5  p        r                        #         @                                                     #FRAGMENT_PLACEMENT%DSIN    #FRAGMENT_PLACEMENT%DCOS    #FRAGMENT_PLACEMENT%TRIM    #FRAGMENT_PLACEMENT%DEXP    #FRAGMENT_PLACEMENT%INT    #FRAGMENT_PLACEMENT%ALLOCATED    #FRAGMENT_PLACEMENT%MAX     #FRAGMENT_PLACEMENT%REAL ¡   #THIS_BOX ¢   #THIS_IM £   #IS ¤   #FRAG_START ¥   #FRAG_TOTAL ¦   #FRAG_ORDER §   #FRAG_PLACED ¨   #THIS_LAMBDA ©   #E_TOTAL ª   #ATTEMPT_PROB «   #NRG_RING_FRAG_TOT ¬   #CBMC_OVERLAP ­   #DEL_OVERLAP ®                                                   DSIN                                                 DCOS                                                 TRIM                                                 DEXP                                                 INT                                                 ALLOCATED                                                  MAX                             @              ¡     REAL           
@ @                               ¢                     
@ @                               £                     
@ @                               ¤                     
                                  ¥                     
                                  ¦                    
                                  §                        p          & p        5 2 r F     5  p        r ¤         5 2 r F     5  p        r ¤                                
D                                 ¨                         p          & p        5 2 r F     5  p        r ¤         5 2 r F     5  p        r ¤                                 
                                 ©     
                
                                ª     
                 
D                                «     
                 
D                                ¬     
                 
D                                 ­                      
D                                 ®            #         @                                   ¯                	   #BUILD_RIGID_FRAGMENT%DEXP °   #BUILD_RIGID_FRAGMENT%INT ±   #BUILD_RIGID_FRAGMENT%ALLOCATED ²   #BUILD_RIGID_FRAGMENT%MAXVAL ³   #BUILD_RIGID_FRAGMENT%MAX ´   #THIS_IM µ   #IS ¶   #THIS_BOX ·   #FRAG_ORDER ¸   #THIS_LAMBDA ¹   #WHICH_ANCHOR º   #ATTEMPT_PROB »   #NRG_RING_FRAG_TOTAL ¼   #CBMC_OVERLAP ½                                              °     DEXP                                            ±     INT                                            ²     ALLOCATED                                            ³     MAXVAL                                            ´     MAX           D @                               µ                      D @                               ¶                      D @                               ·                     D @                               ¸                     
    p          & p        5 2 r F     5  p        r ¶         5 2 r F     5  p        r ¶                                 
  @                              ¹     
                                                 º                      D @                              »     
                 D @                              ¼     
                 
D @                               ½            #         @                                   ¾                   #CUT_REGROW%INT ¿   #CUT_REGROW%ALLOCATED À   #THIS_IM Á   #IS Â   #FRAG_START Ã   #FRAG_END Ä   #FRAG_ORDER Å   #FRAG_TOTAL Æ   #THIS_LAMBDA Ç   #E_PREV È   #ATTEMPT_PROB É   #NRG_RING_FRAG_TOT Ê   #CBMC_OVERLAP Ë   #DEL_OVERLAP Ì                                              ¿     INT                                            À     ALLOCATED           D @                               Á                      D @                               Â                      D @                               Ã                      D @                               Ä                     D @                               Å                         p          & p        5 2 r F     5  p        r Â         5 2 r F     5  p        r Â                                 D @                               Æ                      
  @                              Ç     
                D                                È     
                 D @                              É     
                 D @                              Ê     
                 D @                               Ë                      D @                               Ì            #         @                                  Í                    #IS Î   #FRAG1 Ï   #FRAG2 Ð   #ATOM1 Ñ   #ATOM2 Ò             
                                  Î                     
                                  Ï                     
                                  Ð                     D                                 Ñ                      D                                 Ò            #         @                                  Ó                   #GET_ALIGNER_HANGER%DOT_PRODUCT Ô   #GET_ALIGNER_HANGER%DSQRT Õ   #VEC1 Ö   #VEC2 ×   #ALIGNER Ø   #HANGER Ù                                                                         Ô     DOT_PRODUCT                                            Õ     DSQRT           D @                              Ö                   
 .    p          p            p                                    D @                              ×                   
 /    p          p            p                                    D                                Ø     	              
 2    p          p          p            p          p                                    D                                Ù     	              
 3    p          p          p            p          p                          #         @                                   Ú                   #SINGLE_FRAGMENT_REGROWTH%INT Û   #ALIVE Ü   #IS Ý                                              Û     INT           
@ @                               Ü                     
@ @                               Ý                  ,      fn#fn    Ì   @   j   RUN_VARIABLES       @   J   ENERGY_ROUTINES "   L  @   J   RANDOM_GENERATORS &     @   J   READ_WRITE_CHECKPOINT ,   Ì  Ä       ATOM_CLASS+TYPE_DEFINITIONS 0     H   a   ATOM_CLASS%RXP+TYPE_DEFINITIONS 0   Ø  H   a   ATOM_CLASS%RYP+TYPE_DEFINITIONS 0      H   a   ATOM_CLASS%RZP+TYPE_DEFINITIONS 4   h  H   a   ATOM_CLASS%RXP_NLS+TYPE_DEFINITIONS 4   °  H   a   ATOM_CLASS%RYP_NLS+TYPE_DEFINITIONS 4   ø  H   a   ATOM_CLASS%RZP_NLS+TYPE_DEFINITIONS 4   @  H   a   ATOM_CLASS%RXP_OLD+TYPE_DEFINITIONS 4     H   a   ATOM_CLASS%RYP_OLD+TYPE_DEFINITIONS 4   Ð  H   a   ATOM_CLASS%RZP_OLD+TYPE_DEFINITIONS 2     H   a   ATOM_CLASS%EXIST+TYPE_DEFINITIONS 0   `  `      MOLECULE_CLASS+TYPE_DEFINITIONS >   À  H   a   MOLECULE_CLASS%MOLECULE_TYPE+TYPE_DEFINITIONS 7     H   a   MOLECULE_CLASS%RX_NUM+TYPE_DEFINITIONS 5   P  H   a   MOLECULE_CLASS%LIVE+TYPE_DEFINITIONS 7     H   a   MOLECULE_CLASS%INSIDE+TYPE_DEFINITIONS :   à  H   a   MOLECULE_CLASS%WHICH_BOX+TYPE_DEFINITIONS 5   (  H   a   MOLECULE_CLASS%XCOM+TYPE_DEFINITIONS 5   p  H   a   MOLECULE_CLASS%YCOM+TYPE_DEFINITIONS 5   ¸  H   a   MOLECULE_CLASS%ZCOM+TYPE_DEFINITIONS 7    	  H   a   MOLECULE_CLASS%EULER1+TYPE_DEFINITIONS 7   H	  H   a   MOLECULE_CLASS%EULER2+TYPE_DEFINITIONS 7   	  H   a   MOLECULE_CLASS%EULER3+TYPE_DEFINITIONS 9   Ø	  H   a   MOLECULE_CLASS%XCOM_OLD+TYPE_DEFINITIONS 9    
  H   a   MOLECULE_CLASS%YCOM_OLD+TYPE_DEFINITIONS 9   h
  H   a   MOLECULE_CLASS%ZCOM_OLD+TYPE_DEFINITIONS ;   °
  H   a   MOLECULE_CLASS%EULER1_OLD+TYPE_DEFINITIONS ;   ø
  H   a   MOLECULE_CLASS%EULER2_OLD+TYPE_DEFINITIONS ;   @  H   a   MOLECULE_CLASS%EULER3_OLD+TYPE_DEFINITIONS ;     H   a   MOLECULE_CLASS%CFC_LAMBDA+TYPE_DEFINITIONS 9   Ð  H   a   MOLECULE_CLASS%MAX_DCOM+TYPE_DEFINITIONS =     H   a   MOLECULE_CLASS%MAX_DCOM_OLD+TYPE_DEFINITIONS +   `        BOX_CLASS+TYPE_DEFINITIONS 5   a  P   a   BOX_CLASS%BOX_SHAPE+TYPE_DEFINITIONS 9   ±  H   a   BOX_CLASS%INT_BOX_SHAPE+TYPE_DEFINITIONS 2   ù  ¼   a   BOX_CLASS%LENGTH+TYPE_DEFINITIONS 6   µ  ¼   a   BOX_CLASS%LENGTH_INV+TYPE_DEFINITIONS 5   q  ¼   a   BOX_CLASS%MAX_DELTA+TYPE_DEFINITIONS 3   -  ¼   a   BOX_CLASS%HLENGTH+TYPE_DEFINITIONS 8   é     a   BOX_CLASS%BASIS_LENGTH+TYPE_DEFINITIONS 5        a   BOX_CLASS%COS_ANGLE+TYPE_DEFINITIONS 1   !     a   BOX_CLASS%ANGLE+TYPE_DEFINITIONS 9   ½     a   BOX_CLASS%FACE_DISTANCE+TYPE_DEFINITIONS 2   Y  H   a   BOX_CLASS%VOLUME+TYPE_DEFINITIONS 2   ¡  H   a   BOX_CLASS%DV_MAX+TYPE_DEFINITIONS ,   é  ô       FRAG_CLASS+TYPE_DEFINITIONS 3   Ý  H   a   FRAG_CLASS%NATOMS+TYPE_DEFINITIONS 5   %  H   a   FRAG_CLASS%NCONNECT+TYPE_DEFINITIONS 5   m  H   a   FRAG_CLASS%NANCHORS+TYPE_DEFINITIONS 1   µ  H   a   FRAG_CLASS%TYPE+TYPE_DEFINITIONS 4   ý  H   a   FRAG_CLASS%NCONFIG+TYPE_DEFINITIONS 3   E     a   FRAG_CLASS%ANCHOR+TYPE_DEFINITIONS 2   Ù     a   FRAG_CLASS%ATOMS+TYPE_DEFINITIONS 9   m     a   FRAG_CLASS%FRAG_CONNECT+TYPE_DEFINITIONS 1     H   a   FRAG_CLASS%RING+TYPE_DEFINITIONS 7   I  H   a   FRAG_CLASS%RCUT_VDWSQ+TYPE_DEFINITIONS 8     H   a   FRAG_CLASS%RCUT_COULSQ+TYPE_DEFINITIONS 8   Ù  H   a   FRAG_CLASS%ALPHA_EWALD+TYPE_DEFINITIONS 5   !  n       FRAGMENT_BOND_CLASS+TYPE_DEFINITIONS ?     H   a   FRAGMENT_BOND_CLASS%FRAGMENT1+TYPE_DEFINITIONS ?   ×  H   a   FRAGMENT_BOND_CLASS%FRAGMENT2+TYPE_DEFINITIONS 2     k       FRAG_COORD_CLASS+TYPE_DEFINITIONS 6     H   a   FRAG_COORD_CLASS%RXP+TYPE_DEFINITIONS 6   Ò  H   a   FRAG_COORD_CLASS%RYP+TYPE_DEFINITIONS 6     H   a   FRAG_COORD_CLASS%RZP+TYPE_DEFINITIONS )   b         NFRAGMENTS+RUN_VARIABLES $   î  p       DP+TYPE_DEFINITIONS (   ^  @       KAPPA_INS+RUN_VARIABLES (     @       KAPPA_ROT+RUN_VARIABLES (   Þ  @       KAPPA_DIH+RUN_VARIABLES %            NATOMS+RUN_VARIABLES (   ª  @       CBMC_FLAG+RUN_VARIABLES ,   ê  ¸       MOLECULE_LIST+RUN_VARIABLES (   ¢  Ì       ATOM_LIST+RUN_VARIABLES ,   n  @       GET_FRAGORDER+RUN_VARIABLES (   ®         RRANF+RANDOM_GENERATORS -   /   =      RRANF%IAND+RANDOM_GENERATORS .   l   >      RRANF%ISHFT+RANDOM_GENERATORS -   ª   =      RRANF%IEOR+RANDOM_GENERATORS '   ç   @       DEL_FLAG+RUN_VARIABLES (   '!  ´       FRAG_LIST+RUN_VARIABLES *   Û!  Ò       FRAG_COORDS+RUN_VARIABLES (   ­"  @       IMREPLACE+RUN_VARIABLES (   í"  @       ISREPLACE+RUN_VARIABLES '   -#         BOX_LIST+RUN_VARIABLES (   È#  q       INT_CUBIC+RUN_VARIABLES F   9$  µ       COMPUTE_MOLECULE_NONBOND_INTER_ENERGY+ENERGY_ROUTINES J   î$  <      COMPUTE_MOLECULE_NONBOND_INTER_ENERGY%SUM+ENERGY_ROUTINES I   *%  @   a   COMPUTE_MOLECULE_NONBOND_INTER_ENERGY%IM+ENERGY_ROUTINES I   j%  @   a   COMPUTE_MOLECULE_NONBOND_INTER_ENERGY%IS+ENERGY_ROUTINES R   ª%  @   a   COMPUTE_MOLECULE_NONBOND_INTER_ENERGY%E_INTER_VDW+ENERGY_ROUTINES Q   ê%  @   a   COMPUTE_MOLECULE_NONBOND_INTER_ENERGY%E_INTER_QQ+ENERGY_ROUTINES N   *&  @   a   COMPUTE_MOLECULE_NONBOND_INTER_ENERGY%OVERLAP+ENERGY_ROUTINES #   j&         BETA+RUN_VARIABLES &   ö&  t       MAX_KBT+RUN_VARIABLES =   j'         COMPUTE_RING_FRAGMENT_ENERGY+ENERGY_ROUTINES G   ÷'  @   a   COMPUTE_RING_FRAGMENT_ENERGY%THIS_FRAG+ENERGY_ROUTINES E   7(  @   a   COMPUTE_RING_FRAGMENT_ENERGY%THIS_IM+ENERGY_ROUTINES @   w(  @   a   COMPUTE_RING_FRAGMENT_ENERGY%IS+ENERGY_ROUTINES F   ·(  @   a   COMPUTE_RING_FRAGMENT_ENERGY%THIS_BOX+ENERGY_ROUTINES K   ÷(  @   a   COMPUTE_RING_FRAGMENT_ENERGY%NRG_RING_FRAG+ENERGY_ROUTINES '   7)  ¤       NRG_FRAG+RUN_VARIABLES A   Û)  Ê       COMPUTE_MOLECULE_DIHEDRAL_ENERGY+ENERGY_ROUTINES E   ¥*  <      COMPUTE_MOLECULE_DIHEDRAL_ENERGY%COS+ENERGY_ROUTINES F   á*  =      COMPUTE_MOLECULE_DIHEDRAL_ENERGY%DCOS+ENERGY_ROUTINES J   +  @   a   COMPUTE_MOLECULE_DIHEDRAL_ENERGY%MOLECULE+ENERGY_ROUTINES I   ^+  @   a   COMPUTE_MOLECULE_DIHEDRAL_ENERGY%SPECIES+ENERGY_ROUTINES N   +  @   a   COMPUTE_MOLECULE_DIHEDRAL_ENERGY%ENERGY_DIHED+ENERGY_ROUTINES -   Þ+         FRAGMENT_BONDS+RUN_VARIABLES 1   j,  ½       FRAGMENT_BOND_LIST+RUN_VARIABLES &   '-         ERR_MSG+RUN_VARIABLES $   Ã-  |       TWOPI+RUN_VARIABLES <   ?.        COMPUTE_ATOM_NONBOND_ENERGY+ENERGY_ROUTINES A   U/  =      COMPUTE_ATOM_NONBOND_ENERGY%SQRT+ENERGY_ROUTINES @   /  <      COMPUTE_ATOM_NONBOND_ENERGY%EXP+ENERGY_ROUTINES F   Î/  @   a   COMPUTE_ATOM_NONBOND_ENERGY%THIS_ATOM+ENERGY_ROUTINES J   0  @   a   COMPUTE_ATOM_NONBOND_ENERGY%THIS_MOLECULE+ENERGY_ROUTINES I   N0  @   a   COMPUTE_ATOM_NONBOND_ENERGY%THIS_SPECIES+ENERGY_ROUTINES H   0  @   a   COMPUTE_ATOM_NONBOND_ENERGY%E_INTRA_VDW+ENERGY_ROUTINES H   Î0  @   a   COMPUTE_ATOM_NONBOND_ENERGY%E_INTER_VDW+ENERGY_ROUTINES G   1  @   a   COMPUTE_ATOM_NONBOND_ENERGY%E_INTRA_QQ+ENERGY_ROUTINES G   N1  @   a   COMPUTE_ATOM_NONBOND_ENERGY%E_INTER_QQ+ENERGY_ROUTINES D   1  @   a   COMPUTE_ATOM_NONBOND_ENERGY%OVERLAP+ENERGY_ROUTINES    Î1  ]      BUILD_MOLECULE $   +3  =      BUILD_MOLECULE%DEXP #   h3  <      BUILD_MOLECULE%INT )   ¤3  B      BUILD_MOLECULE%ALLOCATED &   æ3  ?      BUILD_MOLECULE%MAXVAL #   %4  <      BUILD_MOLECULE%MAX '   a4  @   a   BUILD_MOLECULE%THIS_IM "   ¡4  @   a   BUILD_MOLECULE%IS (   á4  @   a   BUILD_MOLECULE%THIS_BOX *   !5  ú   a   BUILD_MOLECULE%FRAG_ORDER +   6  @   a   BUILD_MOLECULE%THIS_LAMBDA ,   [6  @   a   BUILD_MOLECULE%WHICH_ANCHOR ,   6  @   a   BUILD_MOLECULE%ATTEMPT_PROB 3   Û6  @   a   BUILD_MOLECULE%NRG_RING_FRAG_TOTAL ,   7  @   a   BUILD_MOLECULE%CBMC_OVERLAP    [7  Ù       FRAGMENT_ORDER #   48  <      FRAGMENT_ORDER%INT )   p8  B      FRAGMENT_ORDER%ALLOCATED )   ²8  @   a   FRAGMENT_ORDER%THIS_FRAG "   ò8  @   a   FRAGMENT_ORDER%IS *   29  @   a   FRAGMENT_ORDER%FRAG_TOTAL *   r9  ú   a   FRAGMENT_ORDER%FRAG_ORDER $   l:  ú   a   FRAGMENT_ORDER%LIVE '   f;  ú   a   FRAGMENT_ORDER%DEADEND '   `<  ú   a   FRAGMENT_ORDER%CENTRAL #   Z=        FRAGMENT_PLACEMENT (   [?  =      FRAGMENT_PLACEMENT%DSIN (   ?  =      FRAGMENT_PLACEMENT%DCOS (   Õ?  =      FRAGMENT_PLACEMENT%TRIM (   @  =      FRAGMENT_PLACEMENT%DEXP '   O@  <      FRAGMENT_PLACEMENT%INT -   @  B      FRAGMENT_PLACEMENT%ALLOCATED '   Í@  <      FRAGMENT_PLACEMENT%MAX (   	A  =      FRAGMENT_PLACEMENT%REAL ,   FA  @   a   FRAGMENT_PLACEMENT%THIS_BOX +   A  @   a   FRAGMENT_PLACEMENT%THIS_IM &   ÆA  @   a   FRAGMENT_PLACEMENT%IS .   B  @   a   FRAGMENT_PLACEMENT%FRAG_START .   FB  @   a   FRAGMENT_PLACEMENT%FRAG_TOTAL .   B  ú   a   FRAGMENT_PLACEMENT%FRAG_ORDER /   C  ú   a   FRAGMENT_PLACEMENT%FRAG_PLACED /   zD  @   a   FRAGMENT_PLACEMENT%THIS_LAMBDA +   ºD  @   a   FRAGMENT_PLACEMENT%E_TOTAL 0   úD  @   a   FRAGMENT_PLACEMENT%ATTEMPT_PROB 5   :E  @   a   FRAGMENT_PLACEMENT%NRG_RING_FRAG_TOT 0   zE  @   a   FRAGMENT_PLACEMENT%CBMC_OVERLAP /   ºE  @   a   FRAGMENT_PLACEMENT%DEL_OVERLAP %   úE  {      BUILD_RIGID_FRAGMENT *   uG  =      BUILD_RIGID_FRAGMENT%DEXP )   ²G  <      BUILD_RIGID_FRAGMENT%INT /   îG  B      BUILD_RIGID_FRAGMENT%ALLOCATED ,   0H  ?      BUILD_RIGID_FRAGMENT%MAXVAL )   oH  <      BUILD_RIGID_FRAGMENT%MAX -   «H  @   a   BUILD_RIGID_FRAGMENT%THIS_IM (   ëH  @   a   BUILD_RIGID_FRAGMENT%IS .   +I  @   a   BUILD_RIGID_FRAGMENT%THIS_BOX 0   kI  ú   a   BUILD_RIGID_FRAGMENT%FRAG_ORDER 1   eJ  @   a   BUILD_RIGID_FRAGMENT%THIS_LAMBDA 2   ¥J  @   a   BUILD_RIGID_FRAGMENT%WHICH_ANCHOR 2   åJ  @   a   BUILD_RIGID_FRAGMENT%ATTEMPT_PROB 9   %K  @   a   BUILD_RIGID_FRAGMENT%NRG_RING_FRAG_TOTAL 2   eK  @   a   BUILD_RIGID_FRAGMENT%CBMC_OVERLAP    ¥K  2      CUT_REGROW    ×L  <      CUT_REGROW%INT %   M  B      CUT_REGROW%ALLOCATED #   UM  @   a   CUT_REGROW%THIS_IM    M  @   a   CUT_REGROW%IS &   ÕM  @   a   CUT_REGROW%FRAG_START $   N  @   a   CUT_REGROW%FRAG_END &   UN  ú   a   CUT_REGROW%FRAG_ORDER &   OO  @   a   CUT_REGROW%FRAG_TOTAL '   O  @   a   CUT_REGROW%THIS_LAMBDA "   ÏO  @   a   CUT_REGROW%E_PREV (   P  @   a   CUT_REGROW%ATTEMPT_PROB -   OP  @   a   CUT_REGROW%NRG_RING_FRAG_TOT (   P  @   a   CUT_REGROW%CBMC_OVERLAP '   ÏP  @   a   CUT_REGROW%DEL_OVERLAP *   Q  |       GET_COMMON_FRAGMENT_ATOMS -   Q  @   a   GET_COMMON_FRAGMENT_ATOMS%IS 0   ËQ  @   a   GET_COMMON_FRAGMENT_ATOMS%FRAG1 0   R  @   a   GET_COMMON_FRAGMENT_ATOMS%FRAG2 0   KR  @   a   GET_COMMON_FRAGMENT_ATOMS%ATOM1 0   R  @   a   GET_COMMON_FRAGMENT_ATOMS%ATOM2 #   ËR  Ò       GET_ALIGNER_HANGER /   S  D      GET_ALIGNER_HANGER%DOT_PRODUCT )   áS  >      GET_ALIGNER_HANGER%DSQRT (   T     a   GET_ALIGNER_HANGER%VEC1 (   ³T     a   GET_ALIGNER_HANGER%VEC2 +   GU  ´   a   GET_ALIGNER_HANGER%ALIGNER *   ûU  ´   a   GET_ALIGNER_HANGER%HANGER )   ¯V  }       SINGLE_FRAGMENT_REGROWTH -   ,W  <      SINGLE_FRAGMENT_REGROWTH%INT /   hW  @   a   SINGLE_FRAGMENT_REGROWTH%ALIVE ,   ¨W  @   a   SINGLE_FRAGMENT_REGROWTH%IS 