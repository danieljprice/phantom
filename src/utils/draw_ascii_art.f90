!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: draw_ascii_art
!
!  DESCRIPTION:
!  This module contains cool ascii art
!
!  REFERENCES: None
!
!  OWNER: Mark Hutchison
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module draw_ascii_art
 implicit none

 real, parameter :: wait_ms = 20.

 public  :: draw_grainsize_dist_assist,draw_vesta

 private

contains

!---------------------------------------------------------
!--Draws an ascii text: grain-size distribution assistant
!---------------------------------------------------------
subroutine draw_grainsize_dist_assist
 print*,'     @@@@@@@  @@@@@@@   @@@@@@  @@@ @@@  @@@           @@@@@@ @@@ @@@@@@@@ @@@@@@@@      '
 call sleepmilli(wait_ms)
 print*,'    !@@       @@!  @@@ @@!  @@@ @@! @@!@!@@@          !@@     @@!      @@! @@!           '
 call sleepmilli(wait_ms)
 print*,'    !@! @!@!@ @!@!!@!  @!@!@!@! !!@ @!@@!!@! @!@!@!@!  !@@!!  !!@    @!!   @!!!:!        '
 call sleepmilli(wait_ms)
 print*,'    :!!   !!: !!: :!!  !!:  !!! !!: !!:  !!!              !:! !!:  !!:     !!:           '
 call sleepmilli(wait_ms)
 print*,'     :: :: :   :   : :  :   : : :   ::    :           ::.: :  :   :.::.: : : :: ::       '
 call sleepmilli(wait_ms)
 print*,'                                                                                         '
 call sleepmilli(wait_ms)
 print*,'                                                                                         '
 call sleepmilli(wait_ms)
 print*,'@@@@@@@  @@@  @@@@@@ @@@@@@@ @@@@@@@  @@@ @@@@@@@  @@@  @@@ @@@@@@@ @@@  @@@@@@  @@@  @@@'
 call sleepmilli(wait_ms)
 print*,'@@!  @@@ @@! !@@       @!!   @@!  @@@ @@! @@!  @@@ @@!  @@@   @!!   @@! @@!  @@@ @@!@!@@@'
 call sleepmilli(wait_ms)
 print*,'@!@  !@! !!@  !@@!!    @!!   @!@!!@!  !!@ @!@!@!@  @!@  !@!   @!!   !!@ @!@  !@! @!@@!!@!'
 call sleepmilli(wait_ms)
 print*,'!!:  !!! !!:     !:!   !!:   !!: :!!  !!: !!:  !!! !!:  !!!   !!:   !!: !!:  !!! !!:  !!!'
 call sleepmilli(wait_ms)
 print*,':: :  :  :   ::.: :     :     :   : : :   :: : ::   :.:: :     :    :    : :. :  ::    : '
 call sleepmilli(wait_ms)
 print*,'                                                                                         '
 call sleepmilli(wait_ms)
 print*,'                                                                                         '
 call sleepmilli(wait_ms)
 print*,'       @@@@@@   @@@@@@  @@@@@@ @@@  @@@@@@ @@@@@@@  @@@@@@  @@@  @@@ @@@@@@@             '
 call sleepmilli(wait_ms)
 print*,'      @@!  @@@ !@@     !@@     @@! !@@       @!!   @@!  @@@ @@!@!@@@   @!!               '
 call sleepmilli(wait_ms)
 print*,'      @!@!@!@!  !@@!!   !@@!!  !!@  !@@!!    @!!   @!@!@!@! @!@@!!@!   @!!               '
 call sleepmilli(wait_ms)
 print*,'      !!:  !!!     !:!     !:! !!:     !:!   !!:   !!:  !!! !!:  !!!   !!:               '
 call sleepmilli(wait_ms)
 print*,'       :   : : ::.: :  ::.: :  :   ::.: :     :     :   : : ::    :     :                '

end subroutine draw_grainsize_dist_assist


!---------------------------------------------------
!--Draws an ascii picture of Vesta
!---------------------------------------------------
subroutine draw_vesta
 print*,'                                                                                              '
 call sleepmilli(wait_ms)
 print*,'                         ....::        .                                                      '
 call sleepmilli(wait_ms)
 print*,'                      :i,t;,;;.;. ;,,,;,i;i:.                                                 '
 call sleepmilli(wait_ms)
 print*,'                     .:i:i:. .:,:;fti,:ijjtiiji                                               '
 call sleepmilli(wait_ms)
 print*,'               .i;  t:.    ...:i,ii;;ttt;ttttf;tf:                                            '
 call sleepmilli(wait_ms)
 print*,'                 j,,j.;  j  ,;j.:ii,ijf .;it;:ifjjjj                                          '
 call sleepmilli(wait_ms)
 print*,'                .. :fifi,.: .,;i;;jji:;iit;;,,;;ifjffi                                        '
 call sleepmilli(wait_ms)
 print*,'            :      iLLji;   :i;i;;,iii,,,t,;;ii,;,jjjjf;                                      '
 call sleepmilli(wait_ms)
 print*,'            ,:,   ,,Lfjt     t;;i.,,i;tjfjjii;,ii,i;tjtffj                                    '
 call sleepmilli(wait_ms)
 print*,'           :,,,.,;,,iLf     :i;f :,jijtjjtji;i;t: ,;jtjiffjj                                  '
 call sleepmilli(wait_ms)
 print*,'         ,, :   ;i,,:tj     ,;:  ,,tttjititi,t;,;,:jiijfjLfff,.                               '
 call sleepmilli(wait_ms)
 print*,'       :t       ,:;;: t.   :   ., ,i,t;i;ijj.ii;;:jtitjtfLLLLfG                               '
 call sleepmilli(wait_ms)
 print*,'       i         :t,:;,;, :   .  .G, G,i;:tjtt,t:ittftjjLGLGGGGG:                             '
 call sleepmilli(wait_ms)
 print*,'                ,;:i,.,:;;.. ..    .:tifjffjjjt,fjj;ffjjfLLGLDLfLi                            '
 call sleepmilli(wait_ms)
 print*,'               ,::;::j;j .i;..      ,:L;fjjfjffjit: jftjjffLfLDGLfD                           '
 call sleepmilli(wait_ms)
 print*,'             .:.;t.  fi:t,;:: ,.   ;,;;:fLGtLfiLjtfttjjftjGffLffGLfG                          '
 call sleepmilli(wait_ms)
 print*,'         .::: : ,  ,it,ii.t:    i  :i:,jfGjitLtttjfi.fiL;tjjLLLfDGjjG.                        '
 call sleepmilli(wait_ms)
 print*,'         .t   ;    ii;,,tiiit   ,. :;,:LGfji,tjjfLjtLjitt;jjLjLfLGfftEf                       '
 call sleepmilli(wait_ms)
 print*,'        . .   :   , t:,::,ji,i  . .;,;LtGjitGfGDjfGttftjtttiLLGtLLLjGfDL                      '
 call sleepmilli(wait_ms)
 print*,'           ,: i... ,,; ,..t:j,   fL, .DtGfjiGLLGjijfjjjfjjjittLLLGLLGGGDL                     '
 call sleepmilli(wait_ms)
 print*,'          ,..:,:;. :,.. : ,.::;, fj  ,tLGLfti:jjj :GLjtiijjjjfftLLLGLDGDEG                    '
 call sleepmilli(wait_ms)
 print*,'         ,; ,.i.i. ,;.  .. ,:f,   ,ijiLGLLjti;iffjtjjL,itfjjtfjjDfLGGGGDDE:.                  '
 call sleepmilli(wait_ms)
 print*,'          t  ij.;;t, : .  ;j .    t tjDtjGfji,ttj;LjLjjLfLjitffLiEDGfELLGDK.                  '
 call sleepmilli(wait_ms)
 print*,'          .:;: ,,;i,. ,,: ,  .   ,::t,GjfLLii;iijj;tLfjtiftttfffjfGLDLfLGGDE                  '
 call sleepmilli(wait_ms)
 print*,'           ;;  ,,;:::;:i;.. ,   ;,: iLLGGLfjttit;ijfffLtijtLjffftfGDDGfGDGEKE                 '
 call sleepmilli(wait_ms)
 print*,'             :,;i  ii :t  :.;  .;.:,GGGtijLLji,;;iWDfLfjjfGjGfLEjtLtDLLLDDGDKK                '
 call sleepmilli(wait_ms)
 print*,'      :    ..:,:i. L  t,.:i,::,.:.,,LLji iEjLjti;EKDtiDjjjfjf;LfjfjLtGDKGDGGDEK               '
 call sleepmilli(wait_ms)
 print*,'     .:   :.   ,:..f  ;;::,;,,;.:: .tji,:tGLGf;,iGWLt.tjjLLffjffjjGffDDDjLGGGEEG              '
 call sleepmilli(wait_ms)
 print*,'      :   ...::.,;;;;,;,.i::ii,::...,,jitEDGGf,tftEf: fLGGfLLffLfLLLjGGLLfGLDEEEi             '
 call sleepmilli(wait_ms)
 print*,'        . , ,, :;;ii,;;,.,,t;i  i :::;,,iKDfjLjt;jjjffDLGDGGLGfGfGLGfLGfGLDDGEEEK             '
 call sleepmilli(wait_ms)
 print*,'       ..,.,.,ii t::;;:,:,;ti:,t,i :;i,,iWDi.fft;jtjfLDGDGDGLfLLfLtLfLDGGfDLDDEEEK            '
 call sleepmilli(wait_ms)
 print*,'         ; i.:; .t.,jj  :ti:,,jj:: ; ,:,,EGt Lft,tjjijDGLDDLfjffLGffLjLLGfDGLEEEEEE           '
 call sleepmilli(wait_ms)
 print*,'         i.i .,.;; GL    tt;,,;j, ;,,:.j:tfjDLjjtijfjGEGDLGLfjjfKGfjjLLfLLDDGDEEEEK,          '
 call sleepmilli(wait_ms)
 print*,'        : . : ;;::tGL    ititi.;i;.,: iiifjfLftLi:tGiGDDDfGffjjDEGjtjDLGGDGDDGDEEEEK          '
 call sleepmilli(wait_ms)
 print*,'        ::.,  i.::,Di    iiijj:,j.,. .it;Dj;LL;tttjDLDGGLLLLLjfjDGjjiLfEfDGDLEGDDEEKE         '
 call sleepmilli(wait_ms)
 print*,'      .:: :    ;:. Gi   iji;j;.ii,, : tLiLitLj,GfitGDDDLGifffffjLLjijGGLGGGGGDLDLDEKEE        '
 call sleepmilli(wait_ms)
 print*,'     ..:;,.    ;, ..,..tii;;i,ii,; . ;itjjjLfj,Gj;LLDDjLLLtLLLjtjtLLLDGLfGDDGGGLDDDDDKL       '
 call sleepmilli(wait_ms)
 print*,'     :,.ii,    .:..:.,,,ii;;;,ii.,: ;;i;iffffj,jtjfGDDLLGfGGffjLLjLLfGDjjGGDGfDLDDDEEEE:      '
 call sleepmilli(wait_ms)
 print*,'     ;:;t.,    ...::.:;,;;;;,;ii:.::;:jtiLfjjtiftDGfEGjLGfDGLfjLGjjftjLftjGDfLGLLGEDEEKG      '
 call sleepmilli(wait_ms)
 print*,'     ti.,,;     ,. ., i,i:,.,t:;  :f.:ft,jjjjjtifLLDEGfGLDGLfttjLtLjLffffjEGLfLGfGEDEEKE.     '
 call sleepmilli(wait_ms)
 print*,'    :i;::i,.   .,   :.;ji,,;;,:i,.,t j;t;jjjtf,fLfDGGGGGGGGffjtDiLGiLfffLGEGLLLLLDDGEEDKK     '
 call sleepmilli(wait_ms)
 print*,'   ,;,,:.,:;       .:,jji.,.;i; ;.,::;i;fGjttttLjLGGLDGGLLfGttfDtDftffjfEtGELffffLLEEEEKK     '
 call sleepmilli(wait_ms)
 print*,'   ;;,:.:.i    . .:;ti:,it;  t,:,:.t:,jijfj;fjffLLLGGGGGLGLfjtfLGfjifjfLEDGDGGLffGGDDEEWK     '
 call sleepmilli(wait_ms)
 print*,'  ,;,:..j      ,,;iti:  tit  ;;, . ji;,ijfjtEtiijfjDGLGGLLLLLifLfftfjjfjLfLLffffjLDDDEEEKt    '
 call sleepmilli(wait_ms)
 print*,',:, ;, ,;     .:ifji:    ji ,,,,. .f;t,jtGjiLt;:LLLfGLGjfGGffjLLLLfjGttjfLELfjffjGDDDDEEKK    '
 call sleepmilli(wait_ms)
 print*,', ,,i  .      :ttGf       i,;::,.:.;itj;jLj,jjifLfjjtLfLfffLfjLfLLjjttjjfLDGLjLLfLGGDDEKEKi   '
 call sleepmilli(wait_ms)
 print*,'.,,i          ;DGL;       t;,;: :.:,tjjjLDi;tjfLjftjjLGLfjLLLLGGLLtLttjtfLGLffLLLfDLGEEDEEK   '
 call sleepmilli(wait_ms)
 print*,' .,;         .tEEj        ;;;i.:. . tffGtftfijjLLfiGtGjGjGjjLjGLLfjjjijffGGLLLfGLfDLDEDEKKW   '
 call sleepmilli(wait_ms)
 print*,'  ...        ,jDti        i;,;,.:,,:ijjjiLtG;tjijttDGGffLLLGffLLfDGjjtfjjtGDLGLLGGGDEEDKWEKt  '
 call sleepmilli(wait_ms)
 print*,'   .         ,#Et;.       i;;;::,:.iitffjf;ji;jftjjGGGfGfLfGjLffLfGjtjjDttEDLLLGLLGLGGDKWEKW  '
 call sleepmilli(wait_ms)
 print*,'             jKDf.:.      ii;::,,..i;ttjjtttftjjjjGLGLLGtfLjGLjLfjfjLfLfLLDDfLGGLGffGDEKWKWE  '
 call sleepmilli(wait_ms)
 print*,'             Gjjj;:      ,t;i,,,, ,:iijjtjjjfijfjGjtDfGGjjLjfLjGGLjGfGfLtfGGfLLGfGLGLEGKEKWK  '
 call sleepmilli(wait_ms)
 print*,'            . ffft,     :it:,,;:...iftji,jjfjijjG;fjLfGLfjfffLjLffLGjLfjfLLGLDGKDGDGfGGKKWWK  '
 call sleepmilli(wait_ms)
 print*,'               Lft:    .i,,.j;. :.t:Lj,;iiLiijGffijjGEEDtjLfGLLGfffGfjfjjLfDEDGELDDLLGGEKKWW  '
 call sleepmilli(wait_ms)
 print*,'            :   ti  ::it,,;Gf   ..iift,ttfjtijjffifiWEDGttfLLGLffjffffLjfLLGGGGKLfDGLLDDEKW#  '
 call sleepmilli(wait_ms)
 print*,'           ,   ..,.,:itii;:ft   it,tiiiitfLt,LDfttttDKGLt,fGLLGijjjfLLGfLDGGGGGDLLDGGGDEKKWW  '
 call sleepmilli(wait_ms)
 print*,'         .. .:,: : ,::t;t,::t  :t;titf;iitfjfGGtjitifWGfi;GLLDDt;LLffLLLGLjGDDGGDGDGfLDDEKWK  '
 call sleepmilli(wait_ms)
 print*,'         :  ti;,. .:. ,ii,,,,it:iiG:jLi iLf;fGGiittttLDjtfDGfDGi,GLffLfLfGfLGDEDDEDGfGEEKKWL  '
 call sleepmilli(wait_ms)
 print*,'        ...;iii,....:.::,,:.,;i,,iL,jf,:jjj,LjLfiLfjtijLLDLEjfGfiLfjLGGGGDDLGEDDGDGLGLDEKWWf  '
 call sleepmilli(wait_ms)
 print*,'        :;;;,;;,,...:. .,,: ..; tijjjjiififffjjfjfjjjttfGDGDLGGjLLLLGfGDGLGLGDEDGDDGLDDWEWKt  '
 call sleepmilli(wait_ms)
 print*,'      : tt;;,,,:,:.. . ..,;:. :j,ijiftjtL,LLftjjtDGijfftDEELLDjtLGDLLDLEGGGfGDDDDGELLGDKDKWj  '
 call sleepmilli(wait_ms)
 print*,'      . tji.,:,:,:. :   .i;::.t,t,LijiL,jtfjjfjtjDf;t;Gf#LDDLDLLjLGLLLEKfDGLDKGDDGDGLEEEDDWD  '
 call sleepmilli(wait_ms)
 print*,'       jiji :.,::  ,.    ,:;..i::iLitjtijjjjLLfjfLtij;LLDLKDLLEDLfLfDDEKGDELKEGLGDEGEEEEDDKE  '
 call sleepmilli(wait_ms)
 print*,'      ,iiti:.; f.  ;   . .ii: :,;iL,fijjitiLjfGffji;jiGfLDGKDGDEGfjGfDjKEKDLKDfjGDDGKDEDDEKK  '
 call sleepmilli(wait_ms)
 print*,'   .  :ifi;,.;  .  :.,.  ii ,  j.ijjG,tftjt,ttLffj;,fjED#LtELGDGEjfGGjtEDEELEDfjGGGEEEDDDEEK. '
 call sleepmilli(wait_ms)
 print*,'     .ii;t,:,i    . ,,   ;.;: :ji,ttjiit:jt,jiGjLf;:jLEDDLEEGfELDffLLLLELEDEKGfGGLGEDDDDDEEED '
 call sleepmilli(wait_ms)
 print*,'     ;ji:,;.,,    . ,..   :,  ,;;,iifjit;jijtitfffijLGWjGLGGLLEfKGLGfLGDDDEKEDGLGLDEDEGDGDKKK '
 call sleepmilli(wait_ms)
 print*,'     ,ii,i;:;,,., :;.:,  : : i:iittjjjfj;;;f:ijjLjfjEEDLjDLDGLELEGLDfDLKGEEGEDGLLLGEEEGDDEKEK '
 call sleepmilli(wait_ms)
 print*,'     :,iii,,.,,...,,       j .t,,;fGLjLti;;,:fjitfLDEKGftEELGGDLGDfLELfWL#GfKGLfLDEDEEDDDEKKK:'
 call sleepmilli(wait_ms)
 print*,'      ,tj;;:. ,;:::        f :iii,Dt,Lft;i,;,fttifGLDDfjfDGLEGLDDGDL#jGEDWLfEGGLGDEDDDDEDDKKKE'
 call sleepmilli(wait_ms)
 print*,'       ;tii;,i.: :        :tj:t,iiLijjfLi.;,,ijjtjfKfffLLLDLjDGGfDGLDGGEEDDEELDDDDDDDEEKEKEKKK'
 call sleepmilli(wait_ms)
 print*,'       .iitj,;::;         Lij.;;tjtftffLj:,;tLijjjjGfjLfjGDGffLGDEDGWfDEEEEEGDDDELGEEDEEEKKEKW'
 call sleepmilli(wait_ms)
 print*,'        :ttj,.:,:        ;tif  ittj;;jfLf,.tjitjjttjLDLjEGELjfEjDEDDELDEKKEDGDDKGLLEDDEEEKKKWK'
 call sleepmilli(wait_ms)
 print*,'        .tjji, .,        ittL ,;jiDt,ttLf.:jijftttjLLffjDKGfjLDfGGDGDLLEEKLDGDDELLGEEDDEKKDKWK'
 call sleepmilli(wait_ms)
 print*,'        :.ttit,.,       ,iittjj:,fttt;ftj,;ifGL;;jjtfLLjGDGLGGLGDGDGDLLDDKEDDDDKGLDEDEKEEKEKWK'
 call sleepmilli(wait_ms)
 print*,'          itj;;,       it,,iiLL;.j.jtiit,;itGfL,tfj;LGtjKGGfLEfGLGDEGGLGEKDDEEDEGDEEDEEKEEEKKK'
 call sleepmilli(wait_ms)
 print*,'           ;j;;:   :  ,.ij;ijtLi ;tij;;ti;ijjtjttfttLLjjDLGtGEfLEfDDGDDGDEDDDEEEEEEEEKEKKEEKWK'
 call sleepmilli(wait_ms)
 print*,'        :  .;;::  , :t:j,tj:iitj:ji;t;i,ttjftfijtjftfLfffGGfDDjDKfDEGfGDEEEEGEEEDEEEGKEWKKKEWD'
 call sleepmilli(wait_ms)
 print*,'          .  .:...ti ij:iti :jttiti;itti,jijjjjjE.fjjfffLGfDDLEDDDLEGfLGDGGGEEKDDEEEDKKEKEEKWL'
 call sleepmilli(wait_ms)
 print*,'.         : ;   :,:;i;t;t,i jitti;,f .jj:ji:jjjjLjjttffLGGfGDGEDDDLDGGLLDDLDKEEDGEEDEKWEKKEE#L'
 call sleepmilli(wait_ms)
 print*,'         : .:.:  ;i,t,,,;,;,tttjii,L,if;ijtjfttijtiijGfGLGtGDGDDKLGGLfKGGKfDEEWDEDEKEWKEKEEWWD'
 call sleepmilli(wait_ms)
 print*,':   .  .    .j  .;;.;,:;;;it,;iftt;;ttL,ji;fjjiit;iGWGfjtDjGDDDDEGfjffGDDDDDEKEEGEEEKKDKKEEWWf'
 call sleepmilli(wait_ms)
 print*,' .          ::. ,t,i,.. ;i,i,;itii,,tttiftfLjt;;;iKDDLjti;LGEDEDEELLfLLLGDDEEEDDGEKEKDDWEEKWWf'
 call sleepmilli(wait_ms)
 print*,':   ;   .   : L ,i;;:,,; G  ;itfii,,i;tL;;iGft;;t#EGGfjii;tfEEEEEEfLLLGGDGGDEEEDKEKEKDDWEEKW#:'
 call sleepmilli(wait_ms)
 print*,'         :.  :j , i:ii.i j  ijii;iijijt;iiffftiif#DGjjfi;tjiGDWEEDjGGGDLDGEEEEDDKKKEKEEWKKKKW '
 call sleepmilli(wait_ms)
 print*,'         . .:;::;,i;;i: ,t,.,ii;.L.ftjti;jjf,ft;fEGLfjji;tjiDKKDGLLLDDGfGGEKEEGDWDKEWWKWWWKW# '
 call sleepmilli(wait_ms)
 print*,' ..:.    : i. t.t:,.,:; ,j. ij.:,tijtjitftjtjtt;WKGfjLji,:jfDDEfEDfDDDLLjDKEEEDDWEKKD#WWWWW#i '
 call sleepmilli(wait_ms)
 print*,' . : :       ,,:,:i ;,::,,,.:ij;;,,Lij;jLtjLjj,;DKGffjji,,GGGGLfGEDGGGfLjGEDEEEDKKEEEEWWKW##  '
 call sleepmilli(wait_ms)
 print*,'    .:, , ,  :t .ti., ,:.,i.. Ljt,.fitjDijfLLt;jfKGfLfii;fELEfLGDDLGGDLfjLKKDGEEEKEKKWKWWWW#  '
 call sleepmilli(wait_ms)
 print*,'   :  : ::,  ,; ,,i:.:i.:; :;i;tit:tLiiLfttffiffjfGLfjitjEDLEjDDGDGGGLLffLEKGDEKEEKWKWWW#WKi  '
 call sleepmilli(wait_ms)
 print*,'   , ..  ., ;; .,,;ii.,:i:, ;i,f;;j;tGtjfft;ftjLtiLLjjtLGDGGLfLELGEELLfffDDKDKEEKKKWWWWKWW#   '
 call sleepmilli(wait_ms)
 print*,'        ,,,:::: .:;.i:i:;, . L;jjtj;fLjjfjtGfiffj;,t:LfLDLffGLDGKGEEGfjjGDjKEKDEKKWKKWWWWW#   '
 call sleepmilli(wait_ms)
 print*,'      .; j ;, :t.; ;,;.,;:  i,i:jf;tLfjGjjGLf;ttfijfD:jLGLffLGDLKDDDLjjGDGfKEEEKKEWKKWW##W#   '
 call sleepmilli(wait_ms)
 print*,'       ; t..: ,f.  ;t.i.:, .:,i:;D.iLLffijLLjj:fttjGjGLGfffEjGLDEELELffKGGLKEKGEKKKEWW###W;   '
 call sleepmilli(wait_ms)
 print*,'       .ji.,. tf   ;, ii.;:  :;jjG;LGfftLjjGjittjftffLGLtLfGEfLGjEfDLLfDDGLKDWKEEGKWKW####    '
 call sleepmilli(wait_ms)
 print*,'       ;j,.;  jf;   i,:;j :.. j.ftGLLfLjLjffjtiiLGtjLjLjiGGEELDEGDGGfLffDDEE#KKKDDEKKWWW#     '
 call sleepmilli(wait_ms)
 print*,'       ,::;,  jjf   j,.;,    tf.tfLLfLiGLjGjt;;iGLjjLjjjiGDDEDjGDGLELttfLfDDKEKEEDEKK#W#W     '
 call sleepmilli(wait_ms)
 print*,'       .:::t  ;ji;: jt,:.    iD,LiLfLjfLLGGft,.tffGjjjjjffDEKDGfGGGELjtjLLDWWKKEDGDEWWK#      '
 call sleepmilli(wait_ms)
 print*,'        i  i   ,j.,jjf;,.,   jfGijGjffLftffff,::jffGLjtffLtEKKEGGGWGGLfjLfK#W#EWDDGEEW#G      '
 call sleepmilli(wait_ms)
 print*,'        ,  j.  :,.iLjti,;. .t;,LftfjLfGL.fjjEt.:fDjfLjjjLLtDEWLEDLDEGLffLGWKDWKDEDDEK#W       '
 call sleepmilli(wait_ms)
 print*,'         i    : ,ttfiL.,j  :f,;LL,LffLLjtfiiGj;jtjjLtLtLfGjGDKDEEGGEDGGGfEWWKKKDGDEEKW#       '
 call sleepmilli(wait_ms)
 print*,'        .     ;.:;tfiL,;i;.,tiiLf;fffLjjGtitfttjjf,tEf;fjLLGKEDDGDDLLEDGfKKWKKKEDDEEW#.       '
 call sleepmilli(wait_ms)
 print*,'    ,    ,   .: ; ittjiL i,,;ifLtfGLjfffGijjLtj;jj;;DLfffLfGKEGGGEELDEKLjLW#WWKDDDKW##        '
 call sleepmilli(wait_ms)
 print*,'    .    ;.:.   .,,iffiD ;t:i,jitLjjfjLLjtjjLfj:Gf;fLLLffffKEKGDGfDGGLELGK##WEKEEEKW#f        '
 call sleepmilli(wait_ms)
 print*,'        ::G      ;:;,Gjji;i;;t;;tLLtiGGtD;ffftjjL,jjGLELijjDGDLWftLDDLDLLKW#WKKKEKWWW         '
 call sleepmilli(wait_ms)
 print*,'         LD      .i.ti;fti;tiit,,EEf;LLjfjLfjffjLtjfDjtLGfjELGG#GjGDDGLGLD#W#W#LEKKK          '
 call sleepmilli(wait_ms)
 print*,'         fD       i,t;jjf.jL :i.iDDfjftLfGjLftLjftfLfGfGGfGEfGf#LLDEGLDGLGW#EKDDKK#G          '
 call sleepmilli(wait_ms)
 print*,'       .  G       iitj:fj;,j,i.jiEtLjjffjGjjtjLfjtjffLGLLLEDLLLDGDGELLEDGD#WDKWKWKK           '
 call sleepmilli(wait_ms)
 print*,'       .. ,t      j:fL.ii,tfLt,iGijjfjGffGjGitjijfjjKLtLftDGLDLDGDDEGKEDDGWWKDWEKWK           '
 call sleepmilli(wait_ms)
 print*,'         ..t     .iifi;;;i;jLj,iLjLLffffGiiftLjffjLffffEjfKGLGLDKEDEGEEEDDKWWEEWK#            '
 call sleepmilli(wait_ms)
 print*,'     .   .. .    j,iLi,fi:tfGt,ftLjLGLfjfjtjjtLLLtLDiLLtLLKLGGDGEEEEEGDGKGWWWEEE#i            '
 call sleepmilli(wait_ms)
 print*,'           . ,   .L.,;iii,ffLtLfjjLDGjLLLttjjGiLfGLfffD;GDDGjjLEEDLGDDEG#KKWWKKKE             '
 call sleepmilli(wait_ms)
 print*,'            .    i i,;i.::GfDGfjtLGjGGtLjtfitDttjLjLLffjGfDLLfjGKEGLLDKDWDEKW#EW;             '
 call sleepmilli(wait_ms)
 print*,'           ,    ,:  i::f  tLDLittLLtLGLLffjiiLfGLjftEftffDDjfDDLEKKDLKGWEEEGWKKW              '
 call sleepmilli(wait_ms)
 print*,'         .  i   :L  i ;j; :ftGjiiDjfLGGjLjtiijjLGffGLGDLLDjjDDDLGEEKEWKWKKKKWW#               '
 call sleepmilli(wait_ms)
 print*,'            .; , iii  ,i .t.;iti,Kj.LGEf,jtjWL, :fjGGGGfGfLjLffDGLEDjKGKKDEKWW,               '
 call sleepmilli(wait_ms)
 print*,'            .;G ,.jt,, t::ti  ,t;Kt,fWKi;.i#Wjj  fjfLDLLEfffLWLGGLGKK#KEKEK#WL                '
 call sleepmilli(wait_ms)
 print*,'            :,; j ttit;.;,i;  .ttEjf#KEt;.i#KGt:,ftfELtLGGGDLtGDKDGDWWDGDDWDD                 '
 call sleepmilli(wait_ms)
 print*,'            ;;, ;::,tjjt;..,;. ittLDKEDj,ti#ELjiffLfDfLGGLGGGGfGfWDGEGfjEDDK                  '
 call sleepmilli(wait_ms)
 print*,'             i  ,:.,;f,,i  .:i;ti;jKKGDj;fjfELfGLjKLGGGtWfLLGffDLEGEKEGKDKE                   '
 call sleepmilli(wait_ms)
 print*,'               .:.i if..i..i;,t,tfi;GDWLffjDtjtLjjfGfGjfLjGGEEEEDKLWGLjKKK                    '
 call sleepmilli(wait_ms)
 print*,'             :  i  .t,::j  ,,ittii;tLDfGLLLGjjjDjfLDLDLftfLDL;DDLEGGfEDGK                     '
 call sleepmilli(wait_ms)
 print*,'                ,:,f.,,,:,:t,j:;;iitDGLE LfjjjjDjtDDGLGfjLL#E#KKKDDDEfLE                      '
 call sleepmilli(wait_ms)
 print*,'                 i,,.::;, t,;:,;K  jDGKDtjfLjtK:GtLELLEiGELEKGEKEGEGWEE                       '
 call sleepmilli(wait_ms)
 print*,'                .:;i,:.::;:,f.t,D; tDELLfDLDjLGLLffGfLffGEWKWWGEWDDDKG                        '
 call sleepmilli(wait_ms)
 print*,'                 :;:::,ii  ;:,;;;iitDGGLjLLLfGjGGffjjGLGDDDEW#GEKEDKt                         '
 call sleepmilli(wait_ms)
 print*,'                :,:i.:.j;,.:.K  :.;i,EGG;fDDEGLjLGjjGLLGDEW#KDDLWEW                           '
 call sleepmilli(wait_ms)
 print*,'                . ;:  :::  ,,ff ii,;..fD;jjDEEGLDLftLLfDEE#WWDLLGE                            '
 call sleepmilli(wait_ms)
 print*,'                   ..,:,G t. j;:;t ,ttLiGDjffjjGKfLGLLWLiKWtGGLWj                             '
 call sleepmilli(wait_ms)
 print*,'      ;            . ,,,.;,f:tL, t .;DjtfjLjfjtDfLLjL;DGEWLjLfD                               '
 call sleepmilli(wait_ms)
 print*,'       .,i          . .,j,t,,i;;it,,.jfLLffjftLjjGGLLLLK#WfLEG                                '
 call sleepmilli(wait_ms)
 print*,'        ::,           L   ,i,,,;iji :itijtLEGLfffjLfffEEjLGD                                  '
 call sleepmilli(wait_ms)
 print*,'                     , : j:;ittE  i,tji:tj. :jfjLjDLfjfLGf                                    '
 call sleepmilli(wait_ms)
 print*,'           , :        .,:  t ;jtiftLitji;jt jjtjDLLELfji                                      '
 call sleepmilli(wait_ms)
 print*,'          . . :: .   .;    i ,jttLtjtjD   j:jtftffLGLG                                        '
 call sleepmilli(wait_ms)
 print*,'            ;         .     :t,,tLjGfDGj;. .,i.LiLGL                                          '
 call sleepmilli(wait_ms)
 print*,'                         ::  t.,fjLfEDi,j.i,jDDDGL                                            '
 call sleepmilli(wait_ms)
 print*,'                        :i. . ;;,jjt:titjiii;f,                                               '
 call sleepmilli(wait_ms)
 print*,'                       ;,i,;;;,,LtLjtLj;,,:                                                   '
 call sleepmilli(wait_ms)
 print*,'                        i. :;:;jjfii.                                                         '
 call sleepmilli(wait_ms)
 print*,'                       .,;;j;LjjL;                                                            '
 call sleepmilli(wait_ms)
 print*,'                                                                                              '
end subroutine draw_vesta


!---------------------------------------------------
!--Sleep function that accepts millisecond arguments
!---------------------------------------------------
subroutine sleepmilli(dt)
 real, intent(in) :: dt    ! desired sleep interval [ms]
 character(len=100) :: arg ! input argument character string
 integer,dimension(8) :: t ! arguments for date_and_time
 integer :: ms1,ms2        ! start and end times [ms]

 !--Get start time:
 call date_and_time(values=t)
 ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)

 !--Get the command argument, e.g. sleep time in milliseconds:
 call get_command_argument(number=1,value=arg)

 do !--check time:
    call date_and_time(values=t)
    ms2 = (t(5)*3600+t(6)*60+t(7))*1000+t(8)
    if (ms2-ms1 >= dt) exit
 enddo

end subroutine sleepmilli

end module draw_ascii_art
