	  zZ  Õ   k820309              15.0        CÍÜ_                                                                                                           
       io.f90 DECOMP_2D_IO       	       gen@DECOMP_2D_WRITE_ONE gen@DECOMP_2D_READ_ONE gen@DECOMP_2D_WRITE_VAR gen@DECOMP_2D_READ_VAR gen@DECOMP_2D_WRITE_SCALAR gen@DECOMP_2D_READ_SCALAR gen@DECOMP_2D_WRITE_PLANE gen@DECOMP_2D_WRITE_EVERY gen@DECOMP_2D_WRITE_SUBDOMAIN                      @                              
                                                          
                                                              u #WRITE_ONE_REAL    #WRITE_ONE_COMPLEX 
   #MPIIO_WRITE_REAL_COARSE    #         @     @X                                                #WRITE_ONE_REAL%PRESENT    #IPENCIL    #VAR    #FILENAME    #OPT_DECOMP                  @                                PRESENT           
                                                       
@ @                                                 
              &                   &                   &                                                     
  @                                                  1           
 @                                    è             #DECOMP_INFO 	   #         @     @X                             
                   #WRITE_ONE_COMPLEX%PRESENT    #IPENCIL    #VAR    #FILENAME    #OPT_DECOMP                  @                                PRESENT           
                                                       
@ @                                                               &                   &                   &                                                     
  @                                                  1           
 @                                    è             #DECOMP_INFO 	   #         @     @X                                                 #IPENCIL    #VAR    #FILENAME    #ICOARSE              
                                                       
@ @                                                 
 W             &                   &                   &                                                      @                                                   1           
                                                                                                    u #READ_ONE_REAL    #READ_ONE_COMPLEX    #         @     @X                                                #READ_ONE_REAL%PRESENT    #IPENCIL    #VAR    #FILENAME    #OPT_DECOMP                  @                                PRESENT           
                                                       
D @                                                 
 	              &                   &                   &                                                     
  @                                                  1           
 @                                    è             #DECOMP_INFO 	   #         @     @X                                                #READ_ONE_COMPLEX%PRESENT    #IPENCIL    #VAR    #FILENAME    #OPT_DECOMP                   @                                PRESENT           
                                                       
D @                                                                &                   &                   &                                                     
  @                                                  1           
 @                                     è             #DECOMP_INFO 	                                                          u #WRITE_VAR_REAL !   #WRITE_VAR_COMPLEX (   #         @     @X                             !                   #WRITE_VAR_REAL%PRESENT "   #FH #   #DISP $   #IPENCIL %   #VAR &   #OPT_DECOMP '                 @                           "     PRESENT           
@ @                               #                     
D @                              $                      
                                  %                     
@ @                              &                   
              &                   &                   &                                                     
 @                               '     è             #DECOMP_INFO 	   #         @     @X                             (                   #WRITE_VAR_COMPLEX%PRESENT )   #FH *   #DISP +   #IPENCIL ,   #VAR -   #OPT_DECOMP .                 @                           )     PRESENT           
@ @                               *                     
D @                              +                      
                                  ,                     
@ @                              -                                 &                   &                   &                                                     
 @                               .     è             #DECOMP_INFO 	                                                          u #READ_VAR_REAL /   #READ_VAR_COMPLEX 6   #         @     @X                             /                   #READ_VAR_REAL%PRESENT 0   #FH 1   #DISP 2   #IPENCIL 3   #VAR 4   #OPT_DECOMP 5                 @                           0     PRESENT           
@ @                               1                     
D @                              2                      
                                  3                     
D @                              4                   
               &                   &                   &                                                     
 @                               5     è             #DECOMP_INFO 	   #         @     @X                             6                   #READ_VAR_COMPLEX%PRESENT 7   #FH 8   #DISP 9   #IPENCIL :   #VAR ;   #OPT_DECOMP <                 @                           7     PRESENT           
@ @                               8                     
D @                              9                      
                                  :                     
D @                              ;                                  &                   &                   &                                                     
 @                               <     è             #DECOMP_INFO 	                                                          u #WRITE_SCALAR_REAL =   #WRITE_SCALAR_COMPLEX B   #WRITE_SCALAR_INTEGER G   #WRITE_SCALAR_LOGICAL L   #         @     @X                             =                    #FH >   #DISP ?   #N @   #VAR A             
@ @                               >                     
D @                              ?                      
                                  @                    
@ @                              A                    
 !   p          5  p        r @       5  p        r @                     #         @     @X                             B                    #FH C   #DISP D   #N E   #VAR F             
@ @                               C                     
D @                              D                      
                                  E                    
@ @                              F                     "   p          5  p        r E       5  p        r E                     #         @     @X                             G                    #FH H   #DISP I   #N J   #VAR K             
@ @                               H                     
D @                              I                      
                                  J                    
@ @                               K                     #   p          5  p        r J       5  p        r J                     #         @     @X                             L                    #FH M   #DISP N   #N O   #VAR P             
@ @                               M                     
D @                              N                      
                                  O                    
@ @                               P                     $   p          5  p        r O       5  p        r O                                                                            u #READ_SCALAR_REAL Q   #READ_SCALAR_COMPLEX V   #READ_SCALAR_INTEGER [   #READ_SCALAR_LOGICAL `   #         @     @X                             Q                    #FH R   #DISP S   #N T   #VAR U             
@ @                               R                     
D @                              S                      
@ @                               T                    
D @                              U                    
 %    p          5  p        r T       5  p        r T                     #         @     @X                             V                    #FH W   #DISP X   #N Y   #VAR Z             
@ @                               W                     
D @                              X                      
@ @                               Y                    
D @                              Z                     &    p          5  p        r Y       5  p        r Y                     #         @     @X                             [                    #FH \   #DISP ]   #N ^   #VAR _             
@ @                               \                     
D @                              ]                      
@ @                               ^                    
D @                               _                     '    p          5  p        r ^       5  p        r ^                     #         @     @X                             `                    #FH a   #DISP b   #N c   #VAR d             
@ @                               a                     
D @                              b                      
@ @                               c                    
D @                               d                     (    p          5  p        r c       5  p        r c                                                                            u #WRITE_PLANE_3D_REAL e   #WRITE_PLANE_3D_COMPLEX m   #         @     @X                             e                   #WRITE_PLANE_3D_REAL%PRESENT f   #IPENCIL g   #VAR h   #IPLANE i   #N j   #FILENAME k   #OPT_DECOMP l                 @                           f     PRESENT           
                                  g                     
  @                              h                   
 )             &                   &                   &                                                     
                                  i                     
                                  j                     
  @                              k                    1           
 @                               l     è             #DECOMP_INFO 	   #         @     @X                             m                   #WRITE_PLANE_3D_COMPLEX%PRESENT n   #IPENCIL o   #VAR p   #IPLANE q   #N r   #FILENAME s   #OPT_DECOMP t                 @                           n     PRESENT           
                                  o                     
  @                              p                    0             &                   &                   &                                                     
                                  q                     
                                  r                     
  @                              s                    1           
 @                               t     è             #DECOMP_INFO 	                                                          u #WRITE_EVERY_REAL u   #WRITE_EVERY_COMPLEX ~   #         @     @X                             u                   #WRITE_EVERY_REAL%MOD v   #IPENCIL w   #VAR x   #ISKIP y   #JSKIP z   #KSKIP {   #FILENAME |   #FROM1 }                 @                           v     MOD           
                                  w                     
                                 x                   
 7             &                   &                   &                                                     
                                  y                     
                                  z                     
                                  {                     
  @                              |                    1           
                                  }           #         @     @X                             ~                   #WRITE_EVERY_COMPLEX%MOD    #IPENCIL    #VAR    #ISKIP    #JSKIP    #KSKIP    #FILENAME    #FROM1                  @                                MOD           
                                                       
                                                     G             &                   &                   &                                                     
                                                       
                                                       
                                                       
  @                                                  1           
                                                                                                    u #WRITE_SUBDOMAIN    #         @     @X                                              	   #IPENCIL    #VAR    #IS    #IE    #JS    #JE    #KS    #KE    #FILENAME              
                                                       
                                                    
 [             &                   &                   &                                                     
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
  @                                                  1                   @               @           	     'è                   #XST    #XEN    #XSZ    #YST    #YEN    #YSZ    #ZST    #ZEN    #ZSZ    #X1DIST    #Y1DIST    #Y2DIST    #Z2DIST    #X1CNTS    #Y1CNTS    #Y2CNTS     #Z2CNTS ¡   #X1DISP ¢   #Y1DISP £   #Y2DISP ¤   #Z2DISP ¥   #X1COUNT ¦   #Y1COUNT §   #Y2COUNT ¨   #Z2COUNT ©   #EVEN ª                 $                                                              p          p            p                                        $                                                             p          p            p                                        $                                                             p          p            p                                        $                                          $                   p          p            p                                        $                                          0                   p          p            p                                        $                                          <                   p          p            p                                        $                                          H                   p          p            p                                        $                                          T                   p          p            p                                        $                                          `              	     p          p            p                                      $                                          p              
               &                                                       $                                          ¸                             &                                                       $                                                                       &                                                       $                                          H                            &                                                       $                                                                      &                                                       $                                          Ø                            &                                                       $                                                                        &                                                       $                              ¡            h                            &                                                       $                              ¢            °                            &                                                       $                              £            ø                            &                                                       $                              ¤            @                            &                                                       $                              ¥                                        &                                                         $                              ¦     Ð                          $                              §     Ô                          $                              ¨     Ø                          $                              ©     Ü                          $                              ª     à                                                         «                                                                                                      ¬                                                       17#         @                                  ­                    #DECOMP ®                                              ®     è              #DECOMP_INFO 	                                                ¯                                                       22                                            °                                                        ±                                                        ²                         p          p            p                                                                      ³                         p          p            p                                                                      ´                         p          p            p                                                                      µ                         p          p            p                                                                      ¶                         p          p            p                                                                      ·                         p          p            p                                                                      ¸                         p          p            p                                                                      ¹                         p          p            p                                                                      º                         p          p            p                                                                      »                         p          p            p                                                                      ¼                         p          p            p                                                                      ½                         p          p            p                                                                      ¾                         p          p            p                                                                      ¿                         p          p            p                                                                      À                         p          p            p                                                                      Á                         p          p            p                                                                      Â                         p          p            p                                                                      Ã                         p          p            p                                                                      Ä                                                        Å                                                        Æ            #         @                                  Ç                    #ERRORCODE È   #MSG É             
                                  È                     
                                 É                    1                                             Ê                         p          p            p                                       fn#fn "   ¼   õ   b   uapp(DECOMP_2D_IO    ±  @   J  DECOMP_2D    ñ  @   J  MPI (   1         gen@DECOMP_2D_WRITE_ONE    ¹        WRITE_ONE_REAL '   Q  @      WRITE_ONE_REAL%PRESENT '     @   a   WRITE_ONE_REAL%IPENCIL #   Ñ  ¼   a   WRITE_ONE_REAL%VAR (     L   a   WRITE_ONE_REAL%FILENAME *   Ù  Y   a   WRITE_ONE_REAL%OPT_DECOMP "   2        WRITE_ONE_COMPLEX *   Í  @      WRITE_ONE_COMPLEX%PRESENT *     @   a   WRITE_ONE_COMPLEX%IPENCIL &   M  ¼   a   WRITE_ONE_COMPLEX%VAR +   	  L   a   WRITE_ONE_COMPLEX%FILENAME -   U  Y   a   WRITE_ONE_COMPLEX%OPT_DECOMP (   ®  y      MPIIO_WRITE_REAL_COARSE 0   '  @   a   MPIIO_WRITE_REAL_COARSE%IPENCIL ,   g  ¼   a   MPIIO_WRITE_REAL_COARSE%VAR 1   #	  L   a   MPIIO_WRITE_REAL_COARSE%FILENAME 0   o	  @   a   MPIIO_WRITE_REAL_COARSE%ICOARSE '   ¯	  i       gen@DECOMP_2D_READ_ONE    
        READ_ONE_REAL &   ¯
  @      READ_ONE_REAL%PRESENT &   ï
  @   a   READ_ONE_REAL%IPENCIL "   /  ¼   a   READ_ONE_REAL%VAR '   ë  L   a   READ_ONE_REAL%FILENAME )   7  Y   a   READ_ONE_REAL%OPT_DECOMP !           READ_ONE_COMPLEX )   *  @      READ_ONE_COMPLEX%PRESENT )   j  @   a   READ_ONE_COMPLEX%IPENCIL %   ª  ¼   a   READ_ONE_COMPLEX%VAR *   f  L   a   READ_ONE_COMPLEX%FILENAME ,   ²  Y   a   READ_ONE_COMPLEX%OPT_DECOMP (     k       gen@DECOMP_2D_WRITE_VAR    v        WRITE_VAR_REAL '     @      WRITE_VAR_REAL%PRESENT "   R  @   a   WRITE_VAR_REAL%FH $     @   a   WRITE_VAR_REAL%DISP '   Ò  @   a   WRITE_VAR_REAL%IPENCIL #     ¼   a   WRITE_VAR_REAL%VAR *   Î  Y   a   WRITE_VAR_REAL%OPT_DECOMP "   '        WRITE_VAR_COMPLEX *   Æ  @      WRITE_VAR_COMPLEX%PRESENT %     @   a   WRITE_VAR_COMPLEX%FH '   F  @   a   WRITE_VAR_COMPLEX%DISP *     @   a   WRITE_VAR_COMPLEX%IPENCIL &   Æ  ¼   a   WRITE_VAR_COMPLEX%VAR -     Y   a   WRITE_VAR_COMPLEX%OPT_DECOMP '   Û  i       gen@DECOMP_2D_READ_VAR    D        READ_VAR_REAL &   ß  @      READ_VAR_REAL%PRESENT !     @   a   READ_VAR_REAL%FH #   _  @   a   READ_VAR_REAL%DISP &     @   a   READ_VAR_REAL%IPENCIL "   ß  ¼   a   READ_VAR_REAL%VAR )     Y   a   READ_VAR_REAL%OPT_DECOMP !   ô        READ_VAR_COMPLEX )     @      READ_VAR_COMPLEX%PRESENT $   Ò  @   a   READ_VAR_COMPLEX%FH &     @   a   READ_VAR_COMPLEX%DISP )   R  @   a   READ_VAR_COMPLEX%IPENCIL %     ¼   a   READ_VAR_COMPLEX%VAR ,   N  Y   a   READ_VAR_COMPLEX%OPT_DECOMP +   §  ¥       gen@DECOMP_2D_WRITE_SCALAR "   L  j      WRITE_SCALAR_REAL %   ¶  @   a   WRITE_SCALAR_REAL%FH '   ö  @   a   WRITE_SCALAR_REAL%DISP $   6  @   a   WRITE_SCALAR_REAL%N &   v  ´   a   WRITE_SCALAR_REAL%VAR %   *  j      WRITE_SCALAR_COMPLEX (     @   a   WRITE_SCALAR_COMPLEX%FH *   Ô  @   a   WRITE_SCALAR_COMPLEX%DISP '     @   a   WRITE_SCALAR_COMPLEX%N )   T  ´   a   WRITE_SCALAR_COMPLEX%VAR %     j      WRITE_SCALAR_INTEGER (   r  @   a   WRITE_SCALAR_INTEGER%FH *   ²  @   a   WRITE_SCALAR_INTEGER%DISP '   ò  @   a   WRITE_SCALAR_INTEGER%N )   2   ´   a   WRITE_SCALAR_INTEGER%VAR %   æ   j      WRITE_SCALAR_LOGICAL (   P!  @   a   WRITE_SCALAR_LOGICAL%FH *   !  @   a   WRITE_SCALAR_LOGICAL%DISP '   Ð!  @   a   WRITE_SCALAR_LOGICAL%N )   "  ´   a   WRITE_SCALAR_LOGICAL%VAR *   Ä"  ¡       gen@DECOMP_2D_READ_SCALAR !   e#  j      READ_SCALAR_REAL $   Ï#  @   a   READ_SCALAR_REAL%FH &   $  @   a   READ_SCALAR_REAL%DISP #   O$  @   a   READ_SCALAR_REAL%N %   $  ´   a   READ_SCALAR_REAL%VAR $   C%  j      READ_SCALAR_COMPLEX '   ­%  @   a   READ_SCALAR_COMPLEX%FH )   í%  @   a   READ_SCALAR_COMPLEX%DISP &   -&  @   a   READ_SCALAR_COMPLEX%N (   m&  ´   a   READ_SCALAR_COMPLEX%VAR $   !'  j      READ_SCALAR_INTEGER '   '  @   a   READ_SCALAR_INTEGER%FH )   Ë'  @   a   READ_SCALAR_INTEGER%DISP &   (  @   a   READ_SCALAR_INTEGER%N (   K(  ´   a   READ_SCALAR_INTEGER%VAR $   ÿ(  j      READ_SCALAR_LOGICAL '   i)  @   a   READ_SCALAR_LOGICAL%FH )   ©)  @   a   READ_SCALAR_LOGICAL%DISP &   é)  @   a   READ_SCALAR_LOGICAL%N (   )*  ´   a   READ_SCALAR_LOGICAL%VAR *   Ý*  u       gen@DECOMP_2D_WRITE_PLANE $   R+  °      WRITE_PLANE_3D_REAL ,   ,  @      WRITE_PLANE_3D_REAL%PRESENT ,   B,  @   a   WRITE_PLANE_3D_REAL%IPENCIL (   ,  ¼   a   WRITE_PLANE_3D_REAL%VAR +   >-  @   a   WRITE_PLANE_3D_REAL%IPLANE &   ~-  @   a   WRITE_PLANE_3D_REAL%N -   ¾-  L   a   WRITE_PLANE_3D_REAL%FILENAME /   
.  Y   a   WRITE_PLANE_3D_REAL%OPT_DECOMP '   c.  ³      WRITE_PLANE_3D_COMPLEX /   /  @      WRITE_PLANE_3D_COMPLEX%PRESENT /   V/  @   a   WRITE_PLANE_3D_COMPLEX%IPENCIL +   /  ¼   a   WRITE_PLANE_3D_COMPLEX%VAR .   R0  @   a   WRITE_PLANE_3D_COMPLEX%IPLANE )   0  @   a   WRITE_PLANE_3D_COMPLEX%N 0   Ò0  L   a   WRITE_PLANE_3D_COMPLEX%FILENAME 2   1  Y   a   WRITE_PLANE_3D_COMPLEX%OPT_DECOMP *   w1  o       gen@DECOMP_2D_WRITE_EVERY !   æ1  ²      WRITE_EVERY_REAL %   2  <      WRITE_EVERY_REAL%MOD )   Ô2  @   a   WRITE_EVERY_REAL%IPENCIL %   3  ¼   a   WRITE_EVERY_REAL%VAR '   Ð3  @   a   WRITE_EVERY_REAL%ISKIP '   4  @   a   WRITE_EVERY_REAL%JSKIP '   P4  @   a   WRITE_EVERY_REAL%KSKIP *   4  L   a   WRITE_EVERY_REAL%FILENAME '   Ü4  @   a   WRITE_EVERY_REAL%FROM1 $   5  µ      WRITE_EVERY_COMPLEX (   Ñ5  <      WRITE_EVERY_COMPLEX%MOD ,   6  @   a   WRITE_EVERY_COMPLEX%IPENCIL (   M6  ¼   a   WRITE_EVERY_COMPLEX%VAR *   	7  @   a   WRITE_EVERY_COMPLEX%ISKIP *   I7  @   a   WRITE_EVERY_COMPLEX%JSKIP *   7  @   a   WRITE_EVERY_COMPLEX%KSKIP -   É7  L   a   WRITE_EVERY_COMPLEX%FILENAME *   8  @   a   WRITE_EVERY_COMPLEX%FROM1 .   U8  U       gen@DECOMP_2D_WRITE_SUBDOMAIN     ª8        WRITE_SUBDOMAIN (   F9  @   a   WRITE_SUBDOMAIN%IPENCIL $   9  ¼   a   WRITE_SUBDOMAIN%VAR #   B:  @   a   WRITE_SUBDOMAIN%IS #   :  @   a   WRITE_SUBDOMAIN%IE #   Â:  @   a   WRITE_SUBDOMAIN%JS #   ;  @   a   WRITE_SUBDOMAIN%JE #   B;  @   a   WRITE_SUBDOMAIN%KS #   ;  @   a   WRITE_SUBDOMAIN%KE )   Â;  L   a   WRITE_SUBDOMAIN%FILENAME &   <  o      DECOMP_INFO+DECOMP_2D *   }=     a   DECOMP_INFO%XST+DECOMP_2D *   >     a   DECOMP_INFO%XEN+DECOMP_2D *   µ>     a   DECOMP_INFO%XSZ+DECOMP_2D *   Q?     a   DECOMP_INFO%YST+DECOMP_2D *   í?     a   DECOMP_INFO%YEN+DECOMP_2D *   @     a   DECOMP_INFO%YSZ+DECOMP_2D *   %A     a   DECOMP_INFO%ZST+DECOMP_2D *   ÁA     a   DECOMP_INFO%ZEN+DECOMP_2D *   ]B     a   DECOMP_INFO%ZSZ+DECOMP_2D -   ùB     a   DECOMP_INFO%X1DIST+DECOMP_2D -   C     a   DECOMP_INFO%Y1DIST+DECOMP_2D -   !D     a   DECOMP_INFO%Y2DIST+DECOMP_2D -   µD     a   DECOMP_INFO%Z2DIST+DECOMP_2D -   IE     a   DECOMP_INFO%X1CNTS+DECOMP_2D -   ÝE     a   DECOMP_INFO%Y1CNTS+DECOMP_2D -   qF     a   DECOMP_INFO%Y2CNTS+DECOMP_2D -   G     a   DECOMP_INFO%Z2CNTS+DECOMP_2D -   G     a   DECOMP_INFO%X1DISP+DECOMP_2D -   -H     a   DECOMP_INFO%Y1DISP+DECOMP_2D -   ÁH     a   DECOMP_INFO%Y2DISP+DECOMP_2D -   UI     a   DECOMP_INFO%Z2DISP+DECOMP_2D .   éI  H   a   DECOMP_INFO%X1COUNT+DECOMP_2D .   1J  H   a   DECOMP_INFO%Y1COUNT+DECOMP_2D .   yJ  H   a   DECOMP_INFO%Y2COUNT+DECOMP_2D .   ÁJ  H   a   DECOMP_INFO%Z2COUNT+DECOMP_2D +   	K  H   a   DECOMP_INFO%EVEN+DECOMP_2D !   QK  p       MYTYPE+DECOMP_2D $   ÁK  r       REAL_TYPE+DECOMP_2D *   3L  T       GET_DECOMP_INFO+DECOMP_2D 1   L  Y   a   GET_DECOMP_INFO%DECOMP+DECOMP_2D '   àL  r       COMPLEX_TYPE+DECOMP_2D '   RM  @       MYTYPE_BYTES+DECOMP_2D     M  @       NRANK+DECOMP_2D !   ÒM         XSTART+DECOMP_2D    fN         XEND+DECOMP_2D !   úN         YSTART+DECOMP_2D    O         YEND+DECOMP_2D !   "P         ZSTART+DECOMP_2D    ¶P         ZEND+DECOMP_2D    JQ         XSZS+DECOMP_2D    ÞQ         YSZS+DECOMP_2D    rR         ZSZS+DECOMP_2D    S         XSTS+DECOMP_2D    S         YSTS+DECOMP_2D    .T         ZSTS+DECOMP_2D    ÂT         XSZV+DECOMP_2D    VU         YSZV+DECOMP_2D    êU         ZSZV+DECOMP_2D    ~V         XSTV+DECOMP_2D    W         YSTV+DECOMP_2D    ¦W         ZSTV+DECOMP_2D $   :X  @       NX_GLOBAL+DECOMP_2D $   zX  @       NY_GLOBAL+DECOMP_2D $   ºX  @       NZ_GLOBAL+DECOMP_2D *   úX  `       DECOMP_2D_ABORT+DECOMP_2D 4   ZY  @   a   DECOMP_2D_ABORT%ERRORCODE+DECOMP_2D .   Y  L   a   DECOMP_2D_ABORT%MSG+DECOMP_2D     æY         XSIZE+DECOMP_2D 