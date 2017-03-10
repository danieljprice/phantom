%% pasa.cls - version 1.0 
%% Aptara Inc., dated 1st April 2013
%% 
%% For CUP Journals
%%
\NeedsTeXFormat{LaTeX2e}[1997/12/01]
\ProvidesClass{pas}[2013/04/01 v1.0 Publications of the Astronomical Society of Australia (PASA) document class]
%%
\@twosidetrue\@mparswitchtrue\@twocolumnfalse
%%
\DeclareOption{twoside}{\@twosidetrue \@mparswitchtrue}
\DeclareOption{OneCol}{\@twocolumnfalse}
%%
\ExecuteOptions{twoside}
\ProcessOptions
%%
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\bibname{REFERENCES}
\newcommand\indexname{Index}
\newcommand\figurename{Figure}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\chaptername{Chapter}
\newcommand\appendixname{Appendix}
%%
\newcommand\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
%
\newcount\minute
\newcount\hour
\def\currenttime{%
     \minute\time
     \hour\minute
     \divide\hour60
     \the\hour:\multiply\hour60\advance\minute-\hour\the\minute}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boxes & Dimensions %%%%%%%%%%%%%%%%%%%%%
%
\def\jobtag{\@currname}
%
\def\@pushfilename{%
\xdef\@currnamestack{%
	{\@currname}%
	{\@currext}%
	{\the\catcode`\@}%
	\@currnamestack}}%
	\@onlypreamble\@pushfilename%
	\@pushfilename%

\newlength\trimwidth
\newlength\trimheight
\newlength\typewidth
\newlength\typeheight
%
\newdimen\tempdimen
\newdimen\normaltextheight
\newskip\normalbaselineskip
\newdimen\blindfoliodrop
%
\newdimen\enumdim
\newdimen\figheight
\newdimen\figwidth
\newdimen\tabledim
\newskip\tableleftskip
\newskip\tablerightskip
%
\newdimen\mathindent
\newdimen\emathindent
\newdimen\bibindent
%
\newbox\tempbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Boxes & Dimensions %%%%%%%%%%%%%%%%%
%
\def\mdash{\unskip~--\ \ignorespaces}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Headings
\def\rhfont{\fontsize{10}{10}\selectfont\itshape}
\def\foliofont{\fontsize{10}{10}\selectfont}
\def\dropfoliofont{\fontsize{10}{10}\selectfont}
%% Article
\def\sptitlefont{\fontsize{10}{10}\selectfont\bfseries\raggedright\mathversion{bold}}
\def\titlefont{\fontsize{18}{22}\selectfont\raggedright\mathversion{normal}}
\def\authorfont{\fontsize{10}{12}\selectfont\leftskip0pt\rightskip5.5pc plus1fill\mathversion{normal}}
\def\affilfont{\fontsize{8}{9}\selectfont\raggedright\mathversion{normal}}
\def\abstractfont{\fontsize{9}{11}\selectfont\rightskip5.5pc\selectfont}
\def\keywordfont{\fontsize{9}{11}\selectfont\rightskip5.5pc\selectfont}
\def\printdatefont{\fontsize{8}{8}\selectfont}
%% Section
\def\sectionfont{\fontsize{10}{12}\selectfont\bfseries\rightskip1pc plus1fill\mathversion{bold}}
\def\subsectionfont{\fontsize{10}{12}\selectfont\bfseries\rightskip1pc plus1fill\mathversion{bold}}
\def\subsubsectionfont{\fontsize{10}{12}\selectfont\itshape\rightskip1pc plus1fill\mathversion{normal}}
\def\paragraphfont{\fontsize{10}{12}\selectfont\itshape\mathversion{normal}}%
\def\subparagraphfont{\fontsize{10}{12}\selectfont}%
\def\subsubparagraphfont{\fontsize{10}{12}\selectfont}%
%% Figure
\def\figcaptionfont{\fontsize{8}{10}\selectfont\mathversion{normal}}
\def\figcaptionnumfont{\fontsize{8}{10}\selectfont\bfseries}
%% Table
%% BM
\def\indexfont{\footnotesize\raggedright}
\def\bibliofont{\fontsize{9}{11}\selectfont}
%
\def\listfont{}
\def\listdevicefont{}
\def\extractfont{\fontsize{10}{12}\selectfont\leftskip10pt}
%
%%%%%%%%%%%%%%%%%%%%%%%%% Fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\@xvpt{15}
\def\@xviiipt{18}
\def\@xxivpt{24}
%
\normalbaselineskip12pt
%\everydisplay{\fontsize{9}{11}\selectfont}
%
\renewcommand\normalsize{%
   \@setfontsize\normalsize\@xpt{12}
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus.5\p@
   \abovedisplayshortskip \abovedisplayskip
   \belowdisplayshortskip \abovedisplayskip
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}

\def\biggg#1{{\hbox{$\left#1\vbox to20.5\p@{}\right.\n@space$}}}
\def\Biggg#1{{\hbox{$\left#1\vbox to23.5\p@{}\right.\n@space$}}}

\normalsize
%
\newcommand\small{%
   \@setfontsize\small\@ixpt{11}%
   \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
   \setSmallDelims}
%
\def\setSmallDelims{%
\def\big##1{{\hbox{$\left##1\vbox to7.5\p@{}\right.\n@space$}}}%
\def\Big##1{{\hbox{$\left##1\vbox to10.5\p@{}\right.\n@space$}}}%
\def\bigg##1{{\hbox{$\left##1\vbox to13.5\p@{}\right.\n@space$}}}%
\def\Bigg##1{{\hbox{$\left##1\vbox to16.5\p@{}\right.\n@space$}}}%
\def\biggg##1{{\hbox{$\left##1\vbox to19.5\p@{}\right.\n@space$}}}%
\def\Biggg##1{{\hbox{$\left##1\vbox to22.5\p@{}\right.\n@space$}}}%
}
%
\newcommand\footnotesize{%
   \@setfontsize\footnotesize\@viiipt{9}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
   \setFootnotesizeDelims}
%
\def\setFootnotesizeDelims{%
\def\big##1{{\hbox{$\left##1\vbox to6.5\p@{}\right.\n@space$}}}%
\def\Big##1{{\hbox{$\left##1\vbox to9.5\p@{}\right.\n@space$}}}%
\def\bigg##1{{\hbox{$\left##1\vbox to12.5\p@{}\right.\n@space$}}}%
\def\Bigg##1{{\hbox{$\left##1\vbox to15.5\p@{}\right.\n@space$}}}%
\def\biggg##1{{\hbox{$\left##1\vbox to18.5\p@{}\right.\n@space$}}}%
\def\Biggg##1{{\hbox{$\left##1\vbox to21.5\p@{}\right.\n@space$}}}%
}
%   
\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}
\newcommand\Large{\@setfontsize\Large\@xivpt{18}}
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareOldFontCommand{\bi}{\bfseries\itshape}{\bfseries\itshape}
\newcommand{\cal}{\protect\pcal}
\newcommand{\pcal}{\@fontswitch{\relax}{\mathcal}}
\newcommand{\mit}{\protect\pmit}
\newcommand{\pmit}{\@fontswitch{\relax}{\mathnormal}}
%
%%%%%%%%%%%%%%%%% Normal Font %%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\encodingdefault{OT1}%
\fontencoding{OT1}%
%
\DeclareMathSizes{30}{30}{25}{20}
\DeclareMathSizes{18}{18}{15}{11}
\DeclareMathSizes{13}{13}{10}{17}
%
\def\qvec{\mathaccent"0245}
\def\vec#1{\mathaccent"0245{#1}}
\def\lhd{\mathbin{< \hbox to -.43em{}\hbox{\vrule
      \@width .065em \@height .55em \@depth .05em}\hbox to .2em{}}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setlength\trimheight{11in}
\setlength\trimwidth{8.5in}
%
\setlength\typeheight{58pc}
\setlength\typewidth{41.5pc}
%
\setlength\textheight{56pc}
\setlength\textwidth{41.5pc}
%
\setlength\topmargin{.5in}
\setlength\oddsidemargin   {.75in}
\setlength\evensidemargin\trimwidth
\addtolength\evensidemargin{-\textwidth}
\addtolength\evensidemargin{-\oddsidemargin}
%
\setlength\headheight{7.1\p@}% 
\setlength\headsep   {16.5pt}% \typeheight - \textheight - \headheight
\setlength\topskip   {7.4\p@}
%
\newskip\newtextheight
\newtextheight\baselineskip
\advance\newtextheight-\topskip
%
\addtolength\textheight{-\newtextheight}
%
\setlength\normaltextheight{\textheight}

\setlength\footskip{18pt}
\setlength\maxdepth{4\p@}
%
\blindfoliodrop\trimheight
\advance\blindfoliodrop-\typeheight
\advance\blindfoliodrop-\topmargin
\advance\blindfoliodrop-\footskip
\advance\blindfoliodrop18pt
%
\setlength\parindent{10pt}
%
\setlength\marginparwidth {.75in}
\setlength\marginparsep{7\p@}
\setlength\marginparpush{5\p@}
%
\setlength\footnotesep{6.65\p@}
\setlength{\skip\footins}{15.5\p@ \@plus 8\p@ \@minus 8\p@}% change + based on FN font height
\skip\@mpfootins = \skip\footins
%
\setlength\floatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dblfloatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\@fptop{0\p@}
\setlength\@fpsep{8\p@ \@plus 1fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@}
\setlength\@dblfpsep{8\p@ \@plus 1fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
%
\setlength\partopsep{0pt}
\setlength\lineskip{1\p@}% check if it can be flexible
\setlength\normallineskip{1\p@}%
\renewcommand\baselinestretch{}
\setlength\parskip{\z@}%{0\p@ \@plus \p@}
\sloppy
\clubpenalty10000
\widowpenalty10000%
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
%
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Parameters %%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\spreadoutfactor{1}
%
\def\spreadout#1{%
  \gdef\temp{#1}\dimen0 = \spreadoutfactor pt
  \expandafter\dospreadout\temp\endmark\kern-\dimen0}
%
\def\dospreadout{%
  \afterassignment\findospreadout
  \let\next= }
%
\def\findospreadout{%
   \ifx\next\endmark
      \let\nextaction = \relax
   \else
      \let\nextaction = \dospreadout
      \next
      \kern\dimen0
   \fi
   \nextaction}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Article Title  %%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\cleardoublepage{\clearpage\if@twoside \ifodd\c@page\else
    \hbox{}\thispagestyle{empty}\newpage\if@twocolumn\hbox{}\newpage\fi\fi\fi}
%
\def\sptitle#1{\gdef\@sptitle{#1}}
\sptitle{}
%
%%%%%%%%%%%%%%%% New Macro %%%%%%%%%%%
%
\renewcommand\title{\@ifnextchar [{\@titletwo}{\@titleone}}
\def\@titleone#1{\gdef\@title{#1}\gdef\@shorttitle{#1}}
\def\@titletwo[#1]#2{\gdef\@title{\bfseries#2}\gdef\@shorttitle{#1}}
\gdef\@title{\mbox{}}
%%
\def\authormark#1{\xdef\@authormark{#1}}

\def\author{\@ifnextchar [{\@authortwo}{\@authorone}}
\def\@authortwo[#1]#2{\@author{#2}\xdef\@authormark{#1}}
\def\@authorone#1{\@author{#1}\xdef\@authormark{#1}}
%%%
\def\@author#1{\global\advance\authorcount by 1
    \expandafter\gdef\csname author\romannumeral\the\authorcount\endcsname{#1}}
%
\def\affil#1{\expandafter\gdef\csname affil\romannumeral\the\authorcount\endcsname{#1}}
%
\def\and{{\reset@font and}\ }
%
\def\authorsep{\unskip, }
\newcount\authorcount\authorcount=0
\newcount\tempcount
\newcount\authortempcount
%
\def\@authors{%
    \ifnum\authorcount>0%  %%%% Authors available
       \tempcount=1%
       \loop%
       \ifnum\authorcount>0%
          \noindent{\reset@font\authorfont\csname author\romannumeral\the\tempcount \endcsname}\par%
          \noindent{\baselineskip10pt\normalfont\reset@font\affilfont\itshape\csname affil\romannumeral\the\tempcount \endcsname\par}%
           \global\advance\authorcount by -1%
          \advance\tempcount by 1 %
         \vskip7pt
       \repeat
       \global\authortempcount\tempcount
    \fi}%
%%%%%
\newcommand\maketitle{\par%
  \begingroup%
      \twocolumn[\@maketitle]%
%%    \@maketitle
     \renewcommand\thefootnote{\fnsymbol{footnote}}%
     \@thanks
     \endgroup%
  \let\@maketitle\relax%
  \gdef\@authors{}\gdef\@title{}\gdef\thanks{}\gdef\@sptitle{}%
  \renewcommand\thefootnote{\arabic{footnote}}\setcounter{footnote}{0}%
  \@afterheading\@afterindentfalse}
%
\def\@maketitle{%
  \clearpage\null%
  \bgroup
  \setcounter{page}{1}%
    \thispagestyle{plain}%
    \parindent0pt
    \vspace*{5.5pt}
    \hrule height.27pt width 41.5pc\par\vspace*{11.5pt} 
    {\titlefont{\@title\vphantom{pg}}\par}%
    \vskip43.7pt%
    \hrule height.27pt width 41.5pc\par\vspace*{11.3pt} 
    {\@authors\par}\vskip6pt%
    \vskip12pt%
    \box\absbox\par\vskip12pt%
    \box\kwdbox\vspace*{27pt}%
  \egroup}
%
\def\affil#1{{\affilfont #1\vphantom{pqj}\par}}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Article Title  %%%%%%%%%%%%%%%%%%%%%
%
\newbox\absbox
\newenvironment{abstract}{\hsize\textwidth%
  \global\setbox\absbox\vbox\bgroup\par\abstractfont\noindent{\bfseries Abstract}\par\noindent\ignorespaces%
  }{%
  \egroup}
%
\newbox\kwdbox
\newenvironment{keywords}{\hsize\textwidth%
  \global\setbox\kwdbox\vbox\bgroup\par\abstractfont\noindent{\bfseries Keywords:}\ \noindent\ignorespaces%
  }{%
  \egroup}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Sectioning commands  %%%%%%%%%%%%%%%%%%%%
%
\setcounter{secnumdepth}{3}
\newcounter {part}
\newcounter {chapter}
\newcounter {section}[chapter]
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand\thepart          {\Roman{part}}
\renewcommand\thechapter       {\Roman{\c@chapter}}
\renewcommand \thesection {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
%
\def\@seccntformat#1{\csname the#1\endcsname\hskip.5em}
%
\newcommand\section{\@startsection {section}{1}{\z@}{-18pt \@plus -2pt \@minus -.2pt}{6pt}{\sectionfont}}%
\newcommand\subsection{\@startsection{subsection}{2}{\z@}{-18pt \@plus -2pt \@minus -.2pt}{6pt}{\subsectionfont}}%
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}{-10pt \@plus -2pt \@minus -.2pt}{2pt}{\subsubsectionfont}}%
\newcommand\paragraph{\@startsection{paragraph}{4}{\parindent}{-6pt \@plus -2pt \@minus -.2pt}{-.5em}{\paragraphfont}}%
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}{-19pt \@plus -2pt \@minus -.2pt}{-1em}{\subparagraphfont}}%
\newcommand\subsubparagraph{\@startsection{subsubparagraph}{6}{\parindent}{-19pt \@plus -2pt \@minus -.2pt}{-1em}{\subsubparagraphfont}}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Sectioning commands  %%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newskip\topsepi \topsepi12\p@ \@plus2\p@% \@minus.5\p@
\newskip\topsepii \topsepii6pt% \@plus1\p@ 
\newskip\topsepiii \topsepiii6pt% \@plus1\p@ 
\newskip\itemsepi \itemsepi0pt
\newskip\itemsepii \itemsepii0pt
\newskip\itemsepiii \itemsepiii0pt
\newdimen\LabelSep \LabelSep4.5pt
%
\def\@listI{\leftmargin\leftmargini
            \labelwidth\leftmargini
            \advance\labelwidth-\labelsep
            \parsep 0\p@%
            \topsep \topsepi
            \itemsep\itemsepi}%
\let\@listi\@listI
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep\topsepii
              \parsep 0pt
              \itemsep\itemsepii}
\def\@listiii {\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep\topsepiii
              \parsep 0pt
              \itemsep\itemsepiii}
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
%
\setlength\leftmargini  {2.5em}
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv  {1em}
\setlength\leftmarginvi {1em}
\setlength\leftmargin    {\leftmargini}
%
\setlength  \labelsep  {\LabelSep}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
%
\renewcommand\theenumi{\arabic{enumi}}
\renewcommand\theenumii{\alph{enumii}}
\renewcommand\theenumiii{\roman{enumiii}}
\renewcommand\theenumiv{\Alph{enumiv}}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{\theenumii.}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\font\lcir = lcircle10 at 12pt
\newcommand\bulls{\raise1.5pt\hbox{\lcir\char'162}}
\def\textbullet{\leavevmode\raise2pt\hbox{\hskip1pt\bulls}}
\def\textendash{{\bf--}}
\def\textasteriskcentered{\leavevmode\raise-1.5pt\hbox{*}}
\def\textperiodcentered{\leavevmode\raise1.5pt\hbox{\bulls}}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
%
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
%
\newenvironment{verse}
               {\let\\=\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item[]}
               {\endlist}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Quotes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{quote}{\list{}{\rightmargin\leftmargin}%
                \item[]}
               {\endlist}
%
\newenvironment{quotation}{\par\if@nobreak\else\addvspace{12pt plus2pt}\fi%
  \extractfont\parindent0pt\ignorespaces 
}{%
  \par\addvspace{12pt plus2pt}}%
%
\def\leadintextfont{\itshape}
\def\leadintext#1{{\leadintextfont#1}\hskip2.5pt}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Enumerate list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newskip\listtextleftmargin\listtextleftmargin 0pt%24pt
\newskip\listtextleftmarginii\listtextleftmarginii0pt% 24pt
\newskip\listtextleftmarginiii\listtextleftmarginiii0pt% 24pt
\newskip\listtextrightmargin\listtextrightmargin0pt%.5pc
\newskip\listlabelleftskip \listlabelleftskip10pt%3.3pt
\newskip\listlabelleftskipii \listlabelleftskipii0pt%3.3pt
\newskip\listlabelleftskipiii \listlabelleftskipiii0pt%3.3pt
\newskip\abovelistskipi\abovelistskipi12pt plus2pt
\newskip\belowlistskipi\belowlistskipi12pt plus2pt
\newskip\abovelistskipii\abovelistskipii6pt% plus2pt
\newskip\belowlistskipii\belowlistskipii6pt% plus2pt
\newskip\abovelistskipiii\abovelistskipiii6pt% plus2pt
\newskip\belowlistskipiii\belowlistskipiii6pt% plus2pt
\newskip\labelsepi \labelsepi4.5pt
\newskip\labelsepii \labelsepii5.2pt
\newskip\labelsepiii \labelsepiii\z@
\newdimen\enumdimwd
\newif\iflabelrightalign\labelrightaligntrue
%
\def\listdevicefonti{}
\def\listdevicefontii{}
\def\listdevicefontiii{}
\def\listfont{}
\def\textlistlabel{}
%
\def\enummax#1{%
  \labelsep\csname labelsep\romannumeral\the\@enumdepth\endcsname
  \ifdim\listtextleftmargin>\z@\labelsepi0pt\fi
  \ifdim\listtextleftmarginii>\z@\labelsepii0pt\fi
  \ifdim\listtextleftmarginiii>\z@\labelsepiii0pt\fi
  \setbox\tempbox\hbox{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname#1\hskip\labelsep}%
  \enumdim\wd\tempbox
  \setbox\tempbox\hbox{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname#1}%
  \enumdimwd\wd\tempbox
  \expandafter\global\csname leftmargin\romannumeral\the\@enumdepth\endcsname\enumdim
  \ifdim\listtextleftmargin>\z@
     \leftmargini\listtextleftmargin
     \ifdim\listlabelleftskip>\z@
        \advance\leftmargini-\listlabelleftskip
     \fi  
  \fi
  \ifdim\listtextleftmarginii>\z@
     \leftmarginii\listtextleftmarginii
     \ifdim\listlabelleftskipii>\z@
        \advance\leftmarginii-\listlabelleftskipii
     \fi  
  \fi
  \ifdim\listtextleftmarginiii>\z@
     \leftmarginiii\listtextleftmarginiii
     \ifdim\listlabelleftskipiii>\z@
        \advance\leftmarginiii-\listlabelleftskipiii
     \fi  
  \fi
  \ifdim\listlabelleftskip>\z@
     \advance\leftmargini\listlabelleftskip
  \fi
  \ifdim\listlabelleftskipii>\z@
     \advance\leftmarginii\listlabelleftskipii
  \fi
    \ifdim\listlabelleftskipiii>\z@
     \advance\leftmarginiii\listlabelleftskipiii
  \fi}
%
\enummax{1.}
%
\def\enumerate{\@ifnextchar[{\@enumerate}{\@enumerate[\csname label\@enumctr\endcsname]}}%% 
%
\def\@enumerate[#1]{\par
      \ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth\@ne
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \setcounter{\@enumctr}{1}\enummax{#1}%
      \list
          {\csname label\@enumctr\endcsname}{\usecounter{\@enumctr}%
      \topsep\csname abovelistskip\romannumeral\the\@enumdepth\endcsname
      \itemsep\csname itemsep\romannumeral\the\@enumdepth\endcsname
      \listfont %\listparindent18.25pt
      \ifnum \@enumdepth=1 
         \rightmargin\listtextrightmargin
         \advance\rightmargin\rightskip
         \advance\leftmargin\leftskip
         \tempdimen\leftmargini
         \advance\tempdimen-\labelsep
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskip{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi\textlistlabel##1}}}}%
      \else
      \ifnum \@enumdepth=2
         \tempdimen\leftmarginii
         \advance\tempdimen-\labelsep
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskipii{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi##1}}}}%
      \else
         \ifnum \@enumdepth=3
         \tempdimen\leftmarginiii
         \advance\tempdimen-\labelsep
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskipiii{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi##1}}}}%
     \else
         \def\makelabel##1{\hss\llap{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname##1}}%
      \fi\fi
      \fi}\fi}
%
\def\endenumerate{\@topsepadd\csname belowlistskip\romannumeral\the\@enumdepth\endcsname\endlist}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%  unnumlist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{unnumlist}{%
    \list{}{%
%    \listparindent18pt\raggedright
        \leftmargin20pt
        \itemindent-10pt
    \parsep0pt
    \itemsep3pt
    \partopsep0pt}
    \def\makelable##1{##1}%
}{%
    \endlist
}%
%% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% itemize (bullet) %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% define \itemizeiTopSep if it is different from \topsepi
%
\newdimen\itemleftmargin \itemleftmargin11pt % won't be active in enumerate
\newdimen\itemlabelsep \itemlabelsep6.5\p@ % labelsep in itemize for all levels
%
\def\itemize{%
  \ifnum \@itemdepth >3 \@toodeep\else\ifnum\@enumdepth>0\@itemdepth\@enumdepth\fi \advance\@itemdepth \@ne
     \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
     \list{\csname\@itemitem\endcsname}%
        {%
        \topsep\csname abovelistskip\romannumeral\the\@itemdepth\endcsname%\fi
        \itemsep\csname itemsep\romannumeral\the\@itemdepth\endcsname
        \labelsep\itemlabelsep
        \listfont
        \setbox\tempbox\hbox{\csname\@itemitem\endcsname}
        \csname leftmargin\romannumeral\the\@itemdepth\endcsname\wd\tempbox
        \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemlabelsep
        \ifnum\@itemdepth=\@ne\ifnum\@enumdepth=0 % if not inside an enumerate
           \rightmargin\listtextrightmargin
           \advance\rightmargin\rightskip
           \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemleftmargin
           \advance\leftmargini\leftskip
        \else
        \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemleftmargin
           \advance\leftmarginii\leftskip
        \fi\fi
        \leftmargin\csname leftmargin\romannumeral\the\@itemdepth\endcsname
        \def\makelabel##1{\hss\llap{##1}}}
   \fi}%
%
\def\enditemize{%
\@topsepadd\csname belowlistskip\romannumeral\the\@itemdepth\endcsname%\fi
\endlist}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Macro 2: TheoremLike Env. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newif\ifTheorem
\global\Theoremfalse
%
\newcounter{theorem}
\renewcommand\thetheorem{\thechapter.\arabic{theorem}}
%
\newif\ifExample
\global\Examplefalse
%
%
\newif\ifrunon\runonfalse
%
\def\newtheorem#1{\@starthm{#1}\@ifnextchar[{\@othm{#1}}{\@nthm{#1}}}
\def\@starthm#1{\@namedef{#1*}{\@beginstarthm{#1}}%
	\@namedef{end#1*}{\@endstarthm{#1}}}
%
\def\@starthm#1{\@namedef{#1*}{\@beginstarthm{#1}}%
	\@namedef{end#1*}{\@endstarthm}}
%
\def\@beginstarthm#1{\@ifnextchar[{\@opargbeginstarthm{#1}}{\@opargbeginstarthm{#1}[]}}
%
\def\@opargbeginstarthm#1[#2]{\list{}{\topsep12pt plus2pt
 \def\@tempaa{#2}%
	\labelwidth0pt\labelsep9.5pt
	\leftmargin0pt
	\listparindent\parindent
 \edef\@Tempa{\csname #1name\endcsname}
	\def\makelabel##1{##1.}}
      \item[\hskip \labelsep{\csname#1headfont\endcsname \@Tempa{}\ifx\@tempaa\@empty\else\ (#2)\fi}]\csname#1font\endcsname}
%
\def\@endstarthm{\endlist}
%
\def\@xnthm#1#2[#3]{%
  \expandafter\@ifdefinable\csname #1\endcsname
    {\@definecounter{#1}\@newctr{#1}[#3]%
     \expandafter\xdef\csname the#1\endcsname{%
       \expandafter\noexpand\csname the#3\endcsname \@thmcountersep
          \@thmcounter{#1}}%
     \global\@namedef{#1}{\@thm{#1}{#2}{#1}}%
     \global\@namedef{end#1}{\@endtheorem}}}
\def\@ynthm#1#2{%
  \expandafter\@ifdefinable\csname #1\endcsname
    {\@definecounter{#1}%
     \expandafter\xdef\csname the#1\endcsname{\@thmcounter{#1}}%
     \global\@namedef{#1}{\@thm{#1}{#2}{#1}}%
     \global\@namedef{end#1}{\@endtheorem}}}
\def\@othm#1[#2]#3{%
  \@ifundefined{c@#2}{\@nocounterr{#2}}%
    {\expandafter\@ifdefinable\csname #1\endcsname
    {\global\@namedef{the#1}{\@nameuse{the#2}}%
  \global\@namedef{#1}{\@thm{#2}{#3}{#1}}%
  \global\@namedef{end#1}{\@endtheorem}}}}
\def\@thm#1#2#3{%
  \refstepcounter{#1}%
  \@ifnextchar[{\@ythm{#1}{#2}{#3}}{\@xthm{#1}{#2}{#3}}}
%
\def\@xthm#1#2#3{%
  \@opargbegintheorem{#2}{\csname the#1\endcsname}{}{#3}\ignorespaces}
\def\@ythm#1#2#3[#4]{%
  \@opargbegintheorem{#2}{\csname the#1\endcsname}{#4}{#3}\ignorespaces}
%
\def\@endtheorem{\par\addvspace{12pt plus2pt}\endlist}%
%
\def\@opargbegintheorem#1#2#3#4{\par\addvspace{12pt plus2pt}%
%   \renewcommand\labelenumi{{\itshape(\theenumi)}}%
%   \def\labelenumii{{\upshape(\theenumii)}}%
   \def\@tempa{#3}%
   \noindent{\csname #4headfont\endcsname\ifx\@tempa\empty#1\ #2.\hskip9.7pt\else#1\ #2\ {\reset@font\itshape#3}.\hskip9.7pt\fi}%
   \csname #4font\endcsname\ignorespaces}%
%
\def\@endtheorem{\par\addvspace{12pt plus2pt}}%
%
\def\examplename{Example}
\newtheorem{example}{\examplename}[chapter]
\def\theexample{\thesection.\arabic{example}}
\def\examplefont{\itshape}
\def\exampleheadfont{\bfseries}
%
\def\theoremname{Theorem}
\newtheorem{theorem}{\theoremname}[section]
\def\thetheorem{\thesection.\arabic{theorem}}
\def\theoremfont{\itshape}
\def\theoremheadfont{\bfseries}
%
\def\lemmaname{Lemma}
\newtheorem{lemma}{\lemmaname}[section]
\def\thelemma{\thesection.\arabic{lemma}}
\def\lemmafont{\itshape}
\def\lemmaheadfont{\bfseries}
%
\def\propositionname{Proposition}
\newtheorem{proposition}{\propositionname}[section]
\def\theproposition{\thesection.\arabic{proposition}}
\def\propositionfont{\itshape}
\def\propositionheadfont{\bfseries}
%
\def\corollaryname{Corollary}
\newtheorem{corollary}{\corollaryname}[section]
\def\thecorollary{\thesection.\arabic{corollary}}
\def\corollaryfont{\itshape}
\def\corollaryheadfont{\bfseries}
%
\def\definitionname{Def{}inition}
\newtheorem{definition}{\definitionname}[section]
\def\thedefinition{\thesection.\arabic{definition}}
\def\definitionfont{}
\def\definitionheadfont{\bfseries}
%
\def\notationname{Notation}
\newtheorem{notation}{\notationname}[section]
\def\thenotation{\thesection.\arabic{notation}}
\def\notationfont{\itshape}
\def\notationheadfont{\bfseries}
%
\def\assumptionname{Assumption}
\newtheorem{assumption}{\assumptionname}[section]
\def\theassumption{\thesection.\arabic{assumption}}
\def\assumptionfont{}
\def\assumptionheadfont{\bfseries}
%
\def\remarkname{Remark}
\newtheorem{remark}{\remarkname}[section]
\def\theremark{\thesection.\arabic{remark}}
\def\remarkfont{}
\def\remarkheadfont{\bfseries}
%
\newif\ifQEDlogo
\QEDlogofalse
\font\QEDlogofont=msam10 at 9pt
\def\QEDlogo{%
        \ifQEDlogo
 	\else
	  \ifvmode
	    \vskip-\lastskip\hspace*{\fill}%
          \else
  	    \unskip\hskip5pt\ignorespaces
          \fi
          %\par
    	  \hfill{QED}%\QEDlogofont\char'004}\par
          \global\QEDlogofalse
        \fi}
%
\def\prooffont{\reset@font\fontsize{10}{12}\itshape\selectfont\mathversion{normal}}
\font\QEDlogofont=msam10 at 10pt
\def\QEDlogo{\hbox{\QEDlogofont\char'003}}
\def\QEDblogo{\hbox{\QEDlogofont\char'004}}
\newif\ifnologo\nologofalse
\newif\iflogo
\newif\ifblogo\blogofalse
\def\prbox{\unskip\hskip4.5pt\ignorespaces{\fboxsep0pt\fbox{\vrule width5pt height5pt depth0pt}}\global\logofalse}
%\def\QEDlogo{\logotrue{\probox}\global\logofalse}
\newenvironment{proof}{\par%
\leftskip18pt\rightskip\leftskip
			\par%\addvspace{13pt plus2pt}
   \global\QEDlogofalse
   \@ifnextchar[{\@proof}{\@proof[{\prooffont Proof}]}%
   }{%
   \ifblogo\hskip1.2pt
   \QEDblogo
   \else
   \ifnologo
   \else
%   \hfill
%			\QEDlogo
   \fi\fi
   \par%\addvspace{13pt plus2pt}
   }
%
\def\@proof[#1]{\noindent{\prooffont#1\hskip4.8pt}\par\noindent\ignorespaces}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End TheoremLike Env. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Floats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setcounter{topnumber}{5}
\renewcommand\topfraction{.9}
\setcounter{bottomnumber}{5}
\renewcommand\bottomfraction{.9}
\setcounter{totalnumber}{5}
\renewcommand\textfraction{.09}
\renewcommand\floatpagefraction{.901}
\setcounter{dbltopnumber}{5}
\renewcommand\dbltopfraction{.9}
\renewcommand\dblfloatpagefraction{.901}%.901}
%
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{9\p@}
\setlength\belowcaptionskip{5.7\p@}
%
\def\FigName{figure}
%
\long\def\@makecaption#1#2{\figcaptionfont%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{{\bfseries#1}\ #2}%
  \ifdim \wd\@tempboxa >\hsize
    {\bfseries #1}\ #2\par%
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
%%%    \vskip-8pt
  \fi
  \vskip\belowcaptionskip}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newif\ifFP
\global\FPfalse
\newif\ifCO
\global\COfalse
%
\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
%%%
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename\nobreakspace\thefigure.}

\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
%
\def\ArtDir{art/}%
%%\input epsf.sty
%
\newbox\figtempbox
\def\ArtPiece#1{\epsfbox{\ArtDir#1}}%
%
\def\figurebox#1#2{%
    \@ifnextchar[{\@figurebox{#1}{#2}}{\@figurebox{#1}{#2}[]}}
%
\def\@figurebox#1#2[#3]{%
      \gdef\@thirdarg{#3}%
      \ifx\@thirdarg\empty
         \global\figheight=#1
         \global\figwidth=#2
      \else 
         \setbox\figtempbox=\hbox{\epsfbox{\ArtDir\@thirdarg}}%
         \global\figwidth=\wd\figtempbox
         \global\figheight=\ht\figtempbox
      \fi
       \centerline{\figbox}
    \vskip7pt}%
%
\def\figbox{%
     \ifx\@thirdarg\empty
        \vbox{\vskip.3pt\hsize\figwidth
        \hrule\hbox to\figwidth{\vrule\hfill\vbox to\figheight{%
              \hsize\figwidth\vfill}\vrule}\hrule}%
     \else
        \vbox{\hsize\figwidth
 	      \hbox to\figwidth{\vbox to\figheight{\hsize\figwidth\box\figtempbox}}}%
    \fi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcounter{table}[chapter]
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename\space\thetable}%\space

\newenvironment{table}
               {\topsep0pt\@float{table}}
               {\end@float}

\newenvironment{table*}
               {\topsep0pt\@dblfloat{table}}
               {\end@dblfloat}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Floats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Math %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setlength\arraycolsep{1.5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
%
\setlength\columnsep{18\p@}
\setlength\columnseprule{0\p@}
%
\def\thanks#1{\footnotemark
    \protected@xdef\@thanks{\@thanks
        \protect\footnotetext[\the\c@footnote]{#1}}%
}
%
\renewcommand\footnoterule{%
  \kern-5.5\p@
   \hrule width\columnwidth
  \kern 5.5\p@}
\@addtoreset{footnote}{chapter}
\renewcommand\thefootnote{\fnsymbol{footnote}}% \fnsymbol
%
\def\@fnsymbol#1{\ensuremath{\ifcase#1\or *\or \dagger\or \ddagger\or
   \SS\or \|\or \#\or **\or \dagger\dagger\or \ddagger\ddagger 
   \or \SS\SS \or \|\hskip-1pt\|\or \#\#\or ***\or 
   \dagger\dagger\dagger\or \ddagger\ddagger\ddagger\else\@ctrerr\fi\relax}}
%
\newdimen\@footmax
\def\footmax#1{%
  \setbox\tempbox\hbox{\footnotesize$^{#1}$\,}%
  \global\@footmax\wd\tempbox}
%
\footmax{1}
%
\long\def\@makefntext#1{%
  \noindent\@hangfrom{\hbox{\smash{\@makefnmark}}}#1}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Footnote %%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  New Environments %%%%%%%%%%%%%%%%%%%%%

\def\ackfont{\fontsize{9}{11}\selectfont}
\newenvironment{acknowledgements}{%
\section{ACKNOWLEDGEMENTS}
  \noindent\ackfont\ignorespaces%
  }{%
  \par}

\def\ptyfont{\fontsize{10}{12}\rightskip0pt plus1fill\selectfont}
\def\ptysourcefont{\fontsize{10}{12}\leftskip0pt plus1fill\rightskip\poetrytextwidth\selectfont}
%
%\def\extra@space{\vskip\baselineskip}

\newdimen\poetrytextwidth

\newenvironment{poetry}[1][\relax]{\par\addvspace{12pt plus2pt}%26/07
\global\poetrytextwidth\hsize
\global\advance\poetrytextwidth-\wd\tempbox
\let\source\poetrysource
\ptyfont\parindent0pt%
}{%
\par\addvspace{12pt plus2pt}}%
%
\def\poetrysource#1{\par\vskip1.75pt
\ptysourcefont\noindent#1\par}
%
%
\newenvironment{deflist}{\par\addvspace{6pt}%
    \hangindent2em
    \ignorespaces
}{%
\par\addvspace{6pt}}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Page styles  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
\def\doi#1{\gdef\@doi{#1}}
\def\volume#1{\gdef\@volume{#1}}
\def\jnum#1{\gdef\@jnum{#1}}
\def\jslugname#1{\jslugnamefalse\gdef\@jslugname{#1}}
\def\jmonth#1{\gdef\@jmonth{#1}}
\def\jyear#1{\gdef\@jyear{#1}}
\def\jend#1{\gdef\@jend{#1}}
\def\issn#1{\gdef\@issn{#1}}
\def\jid#1{\gdef\@jid{#1}}
\def\jtitle#1{\gdef\@jtitle{#1}}
\def\copyrightstatement#1{\gdef\@copyrightstatement{#1}}
%
\def\jdate#1{\gdef\@jdate{#1}}
%
\def\@volume{}
\def\@jyear{}
\def\@jnum{}
\def\@jslugname{}%
\def\@jmonth{}
\def\@jend{99}
\def\@issn{}
\def\@jid{}
\def\@jtitle{Publications of the Astronomical Society of Australia}
\def\@copyrightstatement{{\copyright}\ Astronomical Society of Australia \the\year; published by Cambridge University Press}
\def\@doi{}
%
\def\endpage#1{\gdef\@endpage{#1}}
%
\def\oddfolio{{\foliofont\thepage}}
\def\evenfolio{{\foliofont\thepage}}
%
\def\ps@headings{%
    \def\@oddfoot{\raggedright\bottomslugline}
    \def\@evenfoot{\raggedright\bottomslugline}
%    
    \def\@evenhead{\evenfolio\hfil{\rhfont\@authormark}}%
    \def\@oddhead{{\rhfont\@shorttitle}\hfil\oddfolio}}
%
\def\ps@empty{%
      \def\@oddfoot{\hfill\lower\blindfoliodrop\hbox{\thepage}\hfill}
      \let\@evenfoot\@oddfoot
      \def\@evenhead{}%
      \def\@oddhead{}%
      \let\@mkboth\markboth
      \let\chaptermark\@gobble
      \let\sectionmark\@gobble%
}
%
\def\sluglinefont{\fontsize{8}{9}\selectfont}

\def\slugline{\vbox{\parindent0pt\sluglinefont\hsize41pc%
    {\itshape \@jtitle}\ifx\@jid\@empty\else\ (\@jid)\fi\par%
     \ifx\@volume\@empty\else,\space Vol. \@volume\par\fi%
     \ifx\@copyrightstatement\@empty\else\@copyrightstatement.\fi\par%
     \ifx\@doi\@empty\else{doi:}\space\@doi.\fi%
\par}
}

\def\bottomslugline{\vtop{\sluglinefont{%
      \ifx\@jid\@empty\else\@jid\ (\@jyear)\fi%
      \ifx\@volume\@empty\else\space\@volume\fi\par}%
      \ifx\@doi\@empty\else{doi:}\@doi\fi}\par}

\def\ps@plain{\let\@mkboth\@gobbletwo%
       \def\@evenfoot{\hfill{\dropfoliofont\thepage}\hfill}%
       \def\@oddfoot{\hfill{\dropfoliofont\thepage}\hfill}%
       \def\@evenhead{\slugline}
       \def\@oddhead{\slugline}}

\def\ps@Empty{%
      \def\@oddfoot{}%
      \let\@evenfoot\@oddfoot
      \def\@evenhead{}%
      \def\@oddhead{}%
      \let\@mkboth\markboth
      \let\chaptermark\@gobble
      \let\sectionmark\@gobble}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Page styles  %%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Bibliography Section   %%%%%%%%%%%%%%%%%%
%
\renewcommand\@biblabel[1]{}% FROM LATEX.LTX
\newcommand\newblock{\hskip .11em \@plus .33em \@minus .07em}

\newcommand\refname{REFERENCES}
\newenvironment{thebibliography}[1]
  {\section*{\refname}%
   \normalfont\small\rmfamily
   \list{}{\usecounter{enumiv}\labelwidth\z@ \leftmargin 9pt \itemindent -9pt}%
   \parindent\z@
   \itemsep0pt\parsep0pt%
   %\parskip -2pt% \@plus .1\p@
   \sloppy\clubpenalty\z@ \widowpenalty\@M
   \sfcode`\.\@m\relax}
  {\def\@noitemerr
   {\@latex@warning{Empty `thebibliography' environment}}%
   \endlist}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Bibliography Section %%%%%%%%%%%%%% 
%% %%
\def\@citex[#1]#2{% FROM LATEX.LTX
  \let\@citea\@empty
  \@cite{\@for\@citeb:=#2\do
    {\@citea\def\@citea{; }%
     \edef\@citeb{\expandafter\@firstofone\@citeb}%
     \if@filesw\immediate\write\@auxout{\string\citation{\@citeb}}\fi
     \@ifundefined{b@\@citeb}{\mbox{\reset@font\bfseries ?}%
       \G@refundefinedtrue
       \@latex@warning
         {Citation `\@citeb' on page \thepage \space undefined}}%
       {\csname b@\@citeb\endcsname}}}{#1}}

%% \def\@citex[#1]#2{%
%%   \let\@citea\@empty
%%   \@cite{\@for\@citeb:=#2\do
%%     {\@citea\def\@citea{; }%\penalty\@m\
%%      \edef\@citeb{\expandafter\@firstofone\@citeb}%\@empty
%%      \if@filesw\immediate\write\@auxout{\string\citation{\@citeb}}\fi
%%      \@ifundefined{b@\@citeb}{\mbox{\reset@font\bfseries ?}%
%%        \G@refundefinedtrue
%%        \@latex@warning
%%          {Citation `\@citeb' on page \thepage \space undefined}}%
%%        {\hbox{\csname b@\@citeb\endcsname}}}}{#1}}


%% %%%\def\@cite#1#2{{\if@tempswa #2\else (#1)\fi}}% FROM LATEX.LTX
\def\@cite#1#2{({#1\if@tempswa , #2\fi})}

\let\@internalcite\cite
\def\cite{\def\citename##1{##1}\@internalcite}

\DeclareRobustCommand\shortcite{\def\citename##1{}\@internalcite}
%% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
\def\@starttoc#1{%
  \begingroup
    \makeatletter
%    \input \jobname.#1%
    \if@filesw
         \expandafter\newwrite\csname tf@#1\endcsname
         \immediate\openout \csname tf@#1\endcsname \jobname.#1\relax
    \fi
    \@nobreakfalse
  \endgroup}
%
 \@starttoc{toc}%
% 

\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg {2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{2}
%
\newcommand\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty{-\@highpenalty}%
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \setlength\@tempdima{3em}%
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hbox to\@pnumwidth{\hss #2}}\par
       \nobreak
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}
    \endgroup
  \fi}

\newcommand\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}

\newcommand\l@section      {\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\l@subsection   {\@dottedtocline{2}{3.8em}{3.2em}}
\newcommand\l@subsubsection{\@dottedtocline{3}{7.0em}{4.1em}}
\newcommand\l@paragraph    {\@dottedtocline{4}{10em}{5em}}
\newcommand\l@subparagraph {\@dottedtocline{5}{12em}{6em}}
\let\l@table\l@figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End FM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Appendix  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\appfont{\fontsize{9}{11}\selectfont}
%
\newcounter{appendix}%
\renewcommand\theappendix{\Alph{appendix}}%
\newenvironment{appendix}{\par\appfont%
  \global\addtocounter{appendix}{1}%
  \setcounter{chapter}{0}%
  %\addtocounter{section}{1}%
  \ifnum\c@appendix>\@ne\else\setcounter{section}{0}\fi%
  \setcounter{table}{0}
  \setcounter{figure}{0}
  \setcounter{equation}{0}%
  \gdef\thechapter{\@Alph{chapter}}
 \renewcommand\thesection{\@Alph\c@section}%
  \renewcommand\thetable{\theappendix\@arabic\c@table}
  \renewcommand\thefigure{\theappendix\@arabic\c@figure}
  \renewcommand\theequation{\theappendix\arabic{equation}}
}{\par}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Index  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\alphabet#1{{\bf#1}\par\vskip.5\baselineskip}
%
\newenvironment{theindex}{\thispagestyle{plain}%
      \columnseprule \z@
      \columnsep 12\p@
      \indexfont\markboth{\indexname}{\indexname}
      \tempdimen\normalbaselineskip
      \advance\tempdimen-\baselineskip
      \twocolumn[\@makeschapterhead{\indexname}\vskip\tempdimen]%
      \parindent\z@
      \let\item\@idxitem
}{%
      \par}
\newcommand\@idxitem  {\par\hangindent3em}
\newcommand\subitem   {\par\hangindent3em\hspace*{1em}}
\newcommand\subsubitem{\par\hangindent3em\hspace*{2em}}
\newcommand\indexspace{\par \vskip\baselineskip}
%
%%%%%%%%%%%%%
%
\def\auand{{and}}
%
\def\tablefont{\fontsize{8}{10}\selectfont}% use plus 1fill if needed
\def\tablecaptionfont{\fontsize{9}{11}\selectfont\leftskip\tableleftskip\rightskip\tablerightskip}% 
\def\tablecaptionnumfont{\fontsize{9}{11}\selectfont\bfseries}% 
\def\TCHfont{\fontsize{8}{9}\selectfont}%
\def\TSCHfont{}%
\def\tabnotefont{\fontsize{8}{9}\selectfont}% use plus 1fill if needed
%
\def\x{@{\extracolsep{\fill}}}
\def\Hline{%
  \noalign{\ifnum0=`}\fi\hrule \@height 1.5pt \futurelet%\@height \arrayrulewidth
   \@tempa\@xhline}
%
\newenvironment{tabnote}{\par\vskip1.5pt% 
  \tabnotefont
}{%
  \par}
%
\def\draftnote{\vbox to 0pt{\vskip-12pt
   \hbox{\ \footnotesize}}%
   \hbox{\ \footnotesize\jobname\vphantom{q}\qquad\jobtag\qquad\today\qquad\currenttime\qquad%
\hfil}}%

\def\Qauthor#1{\marginpar{{\raggedright\footnotesize\bf #1\endgraf}}}
\let\qtoa\Qauthor
%%
%%%%%%%%%%%%%%%%%%%%%%%%% Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%\usepackage{rotating,stfloats-v1}
%
\medmuskip=4mu 
\thickmuskip=5mu 
%
\pagestyle{headings}
\pagenumbering{arabic}      % Arabic page numbers
\sloppy
\frenchspacing
%
\voffset-74.6pt	
\hoffset-72.7pt	
%
\usepackage{graphicx}%,
\usepackage{amsmath,amssymb}
%%\usepackage{txfonts}
%%
\endinput
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Initialization %%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%% End of file `pasa.cls'.
