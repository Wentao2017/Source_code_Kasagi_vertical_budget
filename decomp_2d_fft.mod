	  B#  K   k820309    �          15.0        =~�_                                                                                                           
       fft_generic.f90 DECOMP_2D_FFT              DECOMP_2D_FFT_FINALIZE DECOMP_2D_FFT_GET_SIZE DECOMP_2D_FFT_FORWARD DECOMP_2D_FFT_BACKWARD PHYSICAL_IN_X PHYSICAL_IN_Z gen@DECOMP_2D_FFT_INIT gen@DECOMP_2D_FFT_3D                      @                              
                                                          
                                                              u #FFT_INIT_NOARG    #FFT_INIT_ARG    #FFT_INIT_GENERAL    #         @     @X                                                  #         @     @X                                                #PENCIL              
  @                                          #         @     @X                                                #PENCIL    #NX    #NY 	   #NZ 
             
                                                       
  @                                                    
  @                               	                     
  @                               
                                                                  u #FFT_3D_C2C    #FFT_3D_R2C    #FFT_3D_C2R    #         @     @X                                                 #IN    #OUT    #ISIGN              
                                                                  &                   &                   &                                                     D @                                                                &                   &                   &                                                     
  @                                          #         @     @X                                                 #IN_R    #OUT_C              
  @                                                 
              &                   &                   &                                                     D @                                                                &                   &                   &                                           #         @     @X                                                 #IN_C    #OUT_R              
                                                                  &                   &                   &                                                     D @                                                 
               &                   &                   &                                                             @               @                '�                   #XST    #XEN    #XSZ    #YST    #YEN    #YSZ    #ZST    #ZEN    #ZSZ    #X1DIST    #Y1DIST     #Y2DIST !   #Z2DIST "   #X1CNTS #   #Y1CNTS $   #Y2CNTS %   #Z2CNTS &   #X1DISP '   #Y1DISP (   #Y2DISP )   #Z2DISP *   #X1COUNT +   #Y1COUNT ,   #Y2COUNT -   #Z2COUNT .   #EVEN /                � $                                                              p          p            p                                       � $                                                             p          p            p                                       � $                                                             p          p            p                                       � $                                          $                   p          p            p                                       � $                                          0                   p          p            p                                       � $                                          <                   p          p            p                                       � $                                          H                   p          p            p                                       � $                                          T                   p          p            p                                       � $                                          `              	     p          p            p                                     � $                                          p              
               &                                                      � $                                           �                             &                                                      � $                              !                                         &                                                      � $                              "            H                            &                                                      � $                              #            �                            &                                                      � $                              $            �                            &                                                      � $                              %                                         &                                                      � $                              &            h                            &                                                      � $                              '            �                            &                                                      � $                              (            �                            &                                                      � $                              )            @                            &                                                      � $                              *            �                            &                                                        � $                              +     �                         � $                              ,     �                         � $                              -     �                         � $                              .     �                         � $                              /     �                                                         0                                                         #         @                                  1                    #ERRORCODE 2   #MSG 3             
                                  2                     
                                 3                    1 #         @                                  4                   #DECOMP_INFO_INIT%MOD 5   #DECOMP_INFO_INIT%MAX 6   #DECOMP_INFO_INIT%ALLOCATED 7   #NX 8   #NY 9   #NZ :   #DECOMP ;                 @                           5     MOD               @                           6     MAX               @                           7     ALLOCATED           
                                  8                     
                                  9                     
                                  :                     
                                 ;     �              #DECOMP_INFO    #         @                                  <                    #DECOMP =             
                                 =     �              #DECOMP_INFO                                                >                                                         ?                                          ��������                                                     @                                                      1                                             A                                                      1                                             B                                                      3#         @                                   C                     #         @                                   D                    #ISTART E   #IEND F   #ISIZE G             D                                 E                    	    p          p            p                                    D                                 F                    
    p          p            p                                    D                                 G                        p          p            p                             �   &      fn#fn #   �   �   b   uapp(DECOMP_2D_FFT    y  @   J  DECOMP_2D    �  @   J  GLASSMAN '   �  |       gen@DECOMP_2D_FFT_INIT    u  H      FFT_INIT_NOARG    �  T      FFT_INIT_ARG $     @   a   FFT_INIT_ARG%PENCIL !   Q  l      FFT_INIT_GENERAL (   �  @   a   FFT_INIT_GENERAL%PENCIL $   �  @   a   FFT_INIT_GENERAL%NX $   =  @   a   FFT_INIT_GENERAL%NY $   }  @   a   FFT_INIT_GENERAL%NZ %   �  p       gen@DECOMP_2D_FFT_3D    -  d      FFT_3D_C2C    �  �   a   FFT_3D_C2C%IN    M  �   a   FFT_3D_C2C%OUT !   	  @   a   FFT_3D_C2C%ISIGN    I  ]      FFT_3D_R2C     �  �   a   FFT_3D_R2C%IN_R !   b  �   a   FFT_3D_R2C%OUT_C    	  ]      FFT_3D_C2R     {	  �   a   FFT_3D_C2R%IN_C !   7
  �   a   FFT_3D_C2R%OUT_R &   �
  o      DECOMP_INFO+DECOMP_2D *   b  �   a   DECOMP_INFO%XST+DECOMP_2D *   �  �   a   DECOMP_INFO%XEN+DECOMP_2D *   �  �   a   DECOMP_INFO%XSZ+DECOMP_2D *   6  �   a   DECOMP_INFO%YST+DECOMP_2D *   �  �   a   DECOMP_INFO%YEN+DECOMP_2D *   n  �   a   DECOMP_INFO%YSZ+DECOMP_2D *   
  �   a   DECOMP_INFO%ZST+DECOMP_2D *   �  �   a   DECOMP_INFO%ZEN+DECOMP_2D *   B  �   a   DECOMP_INFO%ZSZ+DECOMP_2D -   �  �   a   DECOMP_INFO%X1DIST+DECOMP_2D -   r  �   a   DECOMP_INFO%Y1DIST+DECOMP_2D -     �   a   DECOMP_INFO%Y2DIST+DECOMP_2D -   �  �   a   DECOMP_INFO%Z2DIST+DECOMP_2D -   .  �   a   DECOMP_INFO%X1CNTS+DECOMP_2D -   �  �   a   DECOMP_INFO%Y1CNTS+DECOMP_2D -   V  �   a   DECOMP_INFO%Y2CNTS+DECOMP_2D -   �  �   a   DECOMP_INFO%Z2CNTS+DECOMP_2D -   ~  �   a   DECOMP_INFO%X1DISP+DECOMP_2D -     �   a   DECOMP_INFO%Y1DISP+DECOMP_2D -   �  �   a   DECOMP_INFO%Y2DISP+DECOMP_2D -   :  �   a   DECOMP_INFO%Z2DISP+DECOMP_2D .   �  H   a   DECOMP_INFO%X1COUNT+DECOMP_2D .     H   a   DECOMP_INFO%Y1COUNT+DECOMP_2D .   ^  H   a   DECOMP_INFO%Y2COUNT+DECOMP_2D .   �  H   a   DECOMP_INFO%Z2COUNT+DECOMP_2D +   �  H   a   DECOMP_INFO%EVEN+DECOMP_2D !   6  p       MYTYPE+DECOMP_2D *   �  `       DECOMP_2D_ABORT+DECOMP_2D 4     @   a   DECOMP_2D_ABORT%ERRORCODE+DECOMP_2D .   F  L   a   DECOMP_2D_ABORT%MSG+DECOMP_2D +   �  �       DECOMP_INFO_INIT+DECOMP_2D 3   R  <      DECOMP_INFO_INIT%MOD+DECOMP_2D=MOD 3   �  <      DECOMP_INFO_INIT%MAX+DECOMP_2D=MAX ?   �  B      DECOMP_INFO_INIT%ALLOCATED+DECOMP_2D=ALLOCATED .     @   a   DECOMP_INFO_INIT%NX+DECOMP_2D .   L  @   a   DECOMP_INFO_INIT%NY+DECOMP_2D .   �  @   a   DECOMP_INFO_INIT%NZ+DECOMP_2D 2   �  Y   a   DECOMP_INFO_INIT%DECOMP+DECOMP_2D /   %  T       DECOMP_INFO_FINALIZE+DECOMP_2D 6   y  Y   a   DECOMP_INFO_FINALIZE%DECOMP+DECOMP_2D     �  @       NRANK+DECOMP_2D &     p       DECOMP_2D_FFT_FORWARD '   �  q       DECOMP_2D_FFT_BACKWARD    �  q       PHYSICAL_IN_X    d   q       PHYSICAL_IN_Z '   �   H       DECOMP_2D_FFT_FINALIZE '   !  i       DECOMP_2D_FFT_GET_SIZE .   �!  �   a   DECOMP_2D_FFT_GET_SIZE%ISTART ,   "  �   a   DECOMP_2D_FFT_GET_SIZE%IEND -   �"  �   a   DECOMP_2D_FFT_GET_SIZE%ISIZE 