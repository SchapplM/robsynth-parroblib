% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:14:42
% EndTime: 2020-08-06 19:14:50
% DurationCPUTime: 8.14s
% Computational Cost: add. (60651->434), mult. (73143->685), div. (8943->6), fcn. (64485->18), ass. (0->296)
t843 = (pkin(2) + pkin(3));
t1009 = -2 * t843;
t842 = (pkin(5) - pkin(6));
t1013 = t842 * t843;
t919 = 2 * pkin(1);
t820 = legFrame(3,3);
t789 = sin(t820);
t792 = cos(t820);
t824 = sin(qJ(1,3));
t830 = cos(qJ(1,3));
t736 = t789 * t824 - t792 * t830;
t823 = sin(qJ(2,3));
t961 = qJ(3,3) * t823;
t773 = pkin(1) + t961;
t785 = t842 * t830;
t743 = t773 * t824 - t785;
t781 = t824 * t842;
t875 = t773 * t830 + t781;
t829 = cos(qJ(2,3));
t931 = t843 * t829;
t714 = t736 * t931 + t743 * t789 - t875 * t792;
t737 = t789 * t830 + t792 * t824;
t715 = t737 * t931 + t743 * t792 + t789 * t875;
t755 = 0.1e1 / (t773 + t931);
t839 = xDP(3);
t840 = xDP(2);
t841 = xDP(1);
t845 = 0.1e1 / qJ(3,3);
t699 = (t823 * t839 + (-t714 * t841 + t715 * t840) * t829 * t755) * t845;
t1012 = -0.2e1 * t699;
t821 = legFrame(2,3);
t790 = sin(t821);
t793 = cos(t821);
t826 = sin(qJ(1,2));
t832 = cos(qJ(1,2));
t738 = t790 * t826 - t793 * t832;
t825 = sin(qJ(2,2));
t962 = qJ(3,2) * t825;
t775 = pkin(1) + t962;
t786 = t842 * t832;
t745 = t775 * t826 - t786;
t782 = t826 * t842;
t873 = t775 * t832 + t782;
t831 = cos(qJ(2,2));
t930 = t843 * t831;
t716 = t738 * t930 + t745 * t790 - t873 * t793;
t739 = t790 * t832 + t793 * t826;
t717 = t739 * t930 + t745 * t793 + t790 * t873;
t756 = 0.1e1 / (t775 + t930);
t847 = 0.1e1 / qJ(3,2);
t700 = (t825 * t839 + (-t716 * t841 + t717 * t840) * t831 * t756) * t847;
t1011 = -0.2e1 * t700;
t822 = legFrame(1,3);
t791 = sin(t822);
t794 = cos(t822);
t828 = sin(qJ(1,1));
t834 = cos(qJ(1,1));
t740 = t791 * t828 - t794 * t834;
t827 = sin(qJ(2,1));
t963 = qJ(3,1) * t827;
t777 = pkin(1) + t963;
t787 = t842 * t834;
t747 = t777 * t828 - t787;
t783 = t828 * t842;
t871 = t777 * t834 + t783;
t833 = cos(qJ(2,1));
t929 = t843 * t833;
t718 = t740 * t929 + t747 * t791 - t871 * t794;
t741 = t791 * t834 + t794 * t828;
t719 = t741 * t929 + t747 * t794 + t791 * t871;
t757 = 0.1e1 / (t777 + t929);
t849 = 0.1e1 / qJ(3,1);
t701 = (t827 * t839 + (-t718 * t841 + t719 * t840) * t833 * t757) * t849;
t1010 = -0.2e1 * t701;
t1008 = (-pkin(5) / 0.2e1 + pkin(6) / 0.2e1) * t1009;
t1007 = t699 * t843;
t1006 = t700 * t843;
t1005 = t701 * t843;
t723 = (-t736 * t840 - t737 * t841) * t755;
t809 = t829 ^ 2;
t1004 = t723 * t809;
t724 = (-t738 * t840 - t739 * t841) * t756;
t811 = t831 ^ 2;
t1003 = t724 * t811;
t725 = (-t740 * t840 - t741 * t841) * t757;
t813 = t833 ^ 2;
t1002 = t725 * t813;
t966 = t823 * pkin(1);
t778 = qJ(3,3) + t966;
t732 = t778 * t824 - t823 * t785;
t772 = pkin(1) + 0.2e1 * t961;
t742 = t772 * t824 - t785;
t866 = t778 * t830 + t823 * t781;
t876 = t772 * t830 + t781;
t941 = (qJ(3,3) + t843) * (-qJ(3,3) + t843);
t886 = t809 * t941;
t708 = -t736 * t886 - (t789 * t742 - t876 * t792) * t931 - (t789 * t732 - t866 * t792) * qJ(3,3);
t709 = t737 * t886 + (t742 * t792 + t876 * t789) * t931 + (t732 * t792 + t866 * t789) * qJ(3,3);
t934 = t843 * t823;
t758 = -qJ(3,3) * t829 + t934;
t945 = t755 * t845;
t687 = t758 * t845 * t839 + (t708 * t841 + t709 * t840) * t945;
t674 = t687 - t1007;
t844 = qJ(3,3) ^ 2;
t666 = t843 * t674 - t844 * t699;
t681 = t699 * t842;
t684 = t687 * t842;
t855 = pkin(1) ^ 2;
t762 = pkin(6) ^ 2 + t855 + (-0.2e1 * pkin(6) + pkin(5)) * pkin(5);
t808 = t829 * t809;
t953 = t723 * t842;
t891 = t823 * t953;
t892 = t723 * t945;
t895 = t699 * t945;
t898 = t687 * t945;
t903 = pkin(1) * qJ(3,3) * t723;
t942 = t762 * t843;
t952 = t723 * t843;
t954 = t723 * t823;
t960 = (t891 - t1007) * t829;
t978 = t809 - 0.1e1;
t982 = -0.3e1 * t844;
t988 = -0.2e1 * t809;
t989 = t843 ^ 2;
t995 = t844 - t989;
t654 = (-(t989 + t982) * t808 * t952 + (-0.3e1 * (-t844 / 0.3e1 + t989) * t961 + t995 * t919) * t1004 + (((-t684 - 0.4e1 * t903) * t843 - t995 * t681) * t823 + t952 * t982 - t723 * t942) * t829) * t892 + ((t666 * t843 + t891 * t941) * t829 + t666 * pkin(1)) * t895 + (pkin(1) * t699 - t960) * t843 * t898 + (((t681 * t1009 + t684) * t809 + (-t762 - t844) * t954 - 0.2e1 * t903 - t684 + t699 * t1008) * t892 + (t666 * t823 + (t988 + 0.1e1) * t723 * t1013) * t895 + (t699 * t934 + t978 * t953) * t898) * qJ(3,3);
t899 = qJ(3,3) * t842 * t809;
t912 = 2 * t843;
t920 = -t842 ^ 2 - t855;
t991 = -2 * pkin(1);
t657 = (-t699 * t899 + (t1013 * t699 - t684) * t823 * t829 + (-t808 * t941 + (t961 * t991 - t844 + t920) * t829 - t773 * t809 * t912) * t723) * t892 + (-t723 * t899 + (t674 + t891) * t931 + t674 * t773) * t895 + (t699 * t773 - t960) * t898;
t985 = 0.2e1 * t823;
t663 = (t758 * t1012 + t687 * t985 + t953) * t723 * t755;
t720 = t723 ^ 2;
t837 = rSges(2,3) + pkin(5);
t835 = pkin(5) + rSges(3,2);
t836 = pkin(2) + rSges(3,1);
t970 = m(3) * t836;
t882 = t835 * t970 - Icges(3,4) - Icges(2,5);
t979 = m(2) * rSges(2,1);
t735 = t837 * t979 + t882;
t838 = m(2) * rSges(2,2);
t927 = Icges(2,6) - Icges(3,6);
t863 = -t837 * t838 + t927;
t817 = -qJ(3,3) - rSges(3,3);
t971 = m(3) * t835;
t908 = t817 * t971;
t726 = -(t863 - t908) * t829 + t823 * t735;
t851 = rSges(2,2) ^ 2;
t853 = rSges(2,1) ^ 2;
t928 = -Icges(2,1) - Icges(3,1);
t868 = Icges(2,2) + Icges(3,3) + (-t851 + t853) * m(2) + t928;
t913 = rSges(3,3) + t836;
t914 = rSges(3,3) - t836;
t729 = -(qJ(3,3) + t913) * (qJ(3,3) + t914) * m(3) + t868;
t883 = -rSges(2,1) * t838 + Icges(2,4) - Icges(3,5);
t752 = -t817 * t970 + t883;
t769 = t970 + t979;
t869 = -(t851 + t853) * m(2) - Icges(3,2) - Icges(2,3);
t850 = rSges(3,3) ^ 2;
t870 = (pkin(2) ^ 2) + t850 + ((2 * pkin(2) + rSges(3,1)) * rSges(3,1));
t974 = m(3) * t817;
t911 = t687 * t974;
t969 = pkin(1) * (t838 + t974);
t990 = 0.2e1 * rSges(3,3);
t998 = qJ(3,3) * t990 + t844;
t923 = t726 * t663 + (-(t870 + t998) * m(3) + t869) * t657 + t654 * t970 + t911 * t1012 + (t752 * t988 + (t729 * t823 + t969) * t829 + t769 * t966 + t752) * t720;
t1001 = t829 * t923;
t965 = t825 * pkin(1);
t779 = qJ(3,2) + t965;
t733 = t779 * t826 - t825 * t786;
t774 = pkin(1) + 0.2e1 * t962;
t744 = t774 * t826 - t786;
t865 = t779 * t832 + t825 * t782;
t874 = t774 * t832 + t782;
t940 = (qJ(3,2) + t843) * (-qJ(3,2) + t843);
t885 = t811 * t940;
t710 = -t738 * t885 - (t790 * t744 - t874 * t793) * t930 - (t790 * t733 - t865 * t793) * qJ(3,2);
t711 = t739 * t885 + (t744 * t793 + t874 * t790) * t930 + (t733 * t793 + t865 * t790) * qJ(3,2);
t933 = t843 * t825;
t759 = -qJ(3,2) * t831 + t933;
t944 = t756 * t847;
t688 = t759 * t847 * t839 + (t710 * t841 + t711 * t840) * t944;
t672 = t688 - t1006;
t846 = qJ(3,2) ^ 2;
t667 = t843 * t672 - t846 * t700;
t682 = t700 * t842;
t685 = t688 * t842;
t810 = t831 * t811;
t950 = t724 * t842;
t889 = t825 * t950;
t890 = t724 * t944;
t894 = t700 * t944;
t897 = t688 * t944;
t904 = pkin(1) * qJ(3,2) * t724;
t949 = t724 * t843;
t951 = t724 * t825;
t959 = (t889 - t1006) * t831;
t977 = t811 - 0.1e1;
t981 = -0.3e1 * t846;
t987 = -0.2e1 * t811;
t994 = t846 - t989;
t655 = (-(t989 + t981) * t810 * t949 + (-0.3e1 * (-t846 / 0.3e1 + t989) * t962 + t994 * t919) * t1003 + (((-t685 - 0.4e1 * t904) * t843 - t994 * t682) * t825 + t949 * t981 - t724 * t942) * t831) * t890 + ((t667 * t843 + t889 * t940) * t831 + t667 * pkin(1)) * t894 + (pkin(1) * t700 - t959) * t843 * t897 + (((t682 * t1009 + t685) * t811 + (-t762 - t846) * t951 - 0.2e1 * t904 - t685 + t700 * t1008) * t890 + (t667 * t825 + (t987 + 0.1e1) * t724 * t1013) * t894 + (t700 * t933 + t977 * t950) * t897) * qJ(3,2);
t900 = qJ(3,2) * t842 * t811;
t658 = (-t700 * t900 + (t1013 * t700 - t685) * t825 * t831 + (-t810 * t940 + (t962 * t991 - t846 + t920) * t831 - t775 * t811 * t912) * t724) * t890 + (-t724 * t900 + (t672 + t889) * t930 + t672 * t775) * t894 + (t700 * t775 - t959) * t897;
t984 = 0.2e1 * t825;
t664 = (t759 * t1011 + t688 * t984 + t950) * t724 * t756;
t721 = t724 ^ 2;
t818 = -qJ(3,2) - rSges(3,3);
t907 = t818 * t971;
t727 = -(t863 - t907) * t831 + t825 * t735;
t730 = -(qJ(3,2) + t913) * (qJ(3,2) + t914) * m(3) + t868;
t753 = -t818 * t970 + t883;
t973 = m(3) * t818;
t910 = t688 * t973;
t968 = pkin(1) * (t838 + t973);
t997 = qJ(3,2) * t990 + t846;
t922 = t727 * t664 + (-(t870 + t997) * m(3) + t869) * t658 + t655 * t970 + t910 * t1011 + (t753 * t987 + (t730 * t825 + t968) * t831 + t769 * t965 + t753) * t721;
t1000 = t831 * t922;
t964 = t827 * pkin(1);
t780 = qJ(3,1) + t964;
t734 = t780 * t828 - t827 * t787;
t776 = pkin(1) + 0.2e1 * t963;
t746 = t776 * t828 - t787;
t864 = t780 * t834 + t827 * t783;
t872 = t776 * t834 + t783;
t939 = (qJ(3,1) + t843) * (-qJ(3,1) + t843);
t884 = t813 * t939;
t712 = -t740 * t884 - (t791 * t746 - t872 * t794) * t929 - (t791 * t734 - t864 * t794) * qJ(3,1);
t713 = t741 * t884 + (t746 * t794 + t872 * t791) * t929 + (t734 * t794 + t864 * t791) * qJ(3,1);
t932 = t843 * t827;
t760 = -qJ(3,1) * t833 + t932;
t943 = t757 * t849;
t689 = t760 * t849 * t839 + (t712 * t841 + t713 * t840) * t943;
t673 = t689 - t1005;
t848 = qJ(3,1) ^ 2;
t668 = t843 * t673 - t848 * t701;
t683 = t701 * t842;
t686 = t689 * t842;
t812 = t833 * t813;
t947 = t725 * t842;
t887 = t827 * t947;
t888 = t725 * t943;
t893 = t701 * t943;
t896 = t689 * t943;
t905 = pkin(1) * qJ(3,1) * t725;
t946 = t725 * t843;
t948 = t725 * t827;
t958 = (t887 - t1005) * t833;
t976 = t813 - 0.1e1;
t980 = -0.3e1 * t848;
t986 = -0.2e1 * t813;
t993 = t848 - t989;
t656 = (-(t989 + t980) * t812 * t946 + (-0.3e1 * (-t848 / 0.3e1 + t989) * t963 + t993 * t919) * t1002 + (((-t686 - 0.4e1 * t905) * t843 - t993 * t683) * t827 + t946 * t980 - t725 * t942) * t833) * t888 + ((t668 * t843 + t887 * t939) * t833 + t668 * pkin(1)) * t893 + (pkin(1) * t701 - t958) * t843 * t896 + (((t683 * t1009 + t686) * t813 + (-t762 - t848) * t948 - 0.2e1 * t905 - t686 + t701 * t1008) * t888 + (t668 * t827 + (t986 + 0.1e1) * t725 * t1013) * t893 + (t701 * t932 + t976 * t947) * t896) * qJ(3,1);
t901 = qJ(3,1) * t842 * t813;
t659 = (-t701 * t901 + (t1013 * t701 - t686) * t827 * t833 + (-t812 * t939 + (t963 * t991 - t848 + t920) * t833 - t777 * t813 * t912) * t725) * t888 + (-t725 * t901 + (t673 + t887) * t929 + t673 * t777) * t893 + (t701 * t777 - t958) * t896;
t983 = 0.2e1 * t827;
t665 = (t760 * t1010 + t689 * t983 + t947) * t725 * t757;
t722 = t725 ^ 2;
t819 = -qJ(3,1) - rSges(3,3);
t906 = t819 * t971;
t728 = -(t863 - t906) * t833 + t827 * t735;
t731 = -(qJ(3,1) + t913) * (qJ(3,1) + t914) * m(3) + t868;
t754 = -t819 * t970 + t883;
t972 = m(3) * t819;
t909 = t689 * t972;
t967 = pkin(1) * (t838 + t972);
t996 = qJ(3,1) * t990 + t848;
t921 = t728 * t665 + (-(t870 + t996) * m(3) + t869) * t659 + t656 * t970 + t909 * t1010 + (t754 * t986 + (t731 * t827 + t967) * t833 + t769 * t964 + t754) * t722;
t999 = t833 * t921;
t975 = m(2) * t837;
t938 = t823 * t835;
t937 = t825 * t835;
t936 = t827 * t835;
t693 = t752 * t699;
t696 = t699 ^ 2;
t751 = rSges(2,1) * t975 + t882;
t761 = t769 * t919;
t862 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - (t837 ^ 2 + t851 + t855) * m(2) - Icges(1,3) + t928;
t880 = t835 ^ 2 + t850 + t855;
t881 = rSges(2,2) * t975 - t927;
t917 = 0.2e1 * t969;
t918 = -0.2e1 * t971;
t926 = -(-t729 * t809 - (t752 * t985 + t761) * t829 + t823 * t917 - (t880 + t998) * m(3) + t862) * t663 - t726 * t657 + m(3) * t654 * t938 + 0.4e1 * (-t693 - t911 / 0.2e1) * t1004 - (-0.2e1 * (-t687 * t970 + t699 * t729) * t954 - t699 * (t687 * t918 + t699 * t751 + t723 * t917)) * t829 - ((t881 + t908) * t696 + (m(3) * t687 - t699 * t769) * t723 * t919) * t823 - 0.2e1 * (-t693 - t911) * t723;
t694 = t753 * t700;
t697 = t700 ^ 2;
t916 = 0.2e1 * t968;
t925 = -(-t730 * t811 - (t753 * t984 + t761) * t831 + t825 * t916 - (t880 + t997) * m(3) + t862) * t664 - t727 * t658 + m(3) * t655 * t937 + 0.4e1 * (-t694 - t910 / 0.2e1) * t1003 - (-0.2e1 * (-t688 * t970 + t700 * t730) * t951 - t700 * (t688 * t918 + t700 * t751 + t724 * t916)) * t831 - ((t881 + t907) * t697 + (m(3) * t688 - t700 * t769) * t724 * t919) * t825 - 0.2e1 * (-t694 - t910) * t724;
t695 = t754 * t701;
t698 = t701 ^ 2;
t915 = 0.2e1 * t967;
t924 = -(-t731 * t813 - (t754 * t983 + t761) * t833 + t827 * t915 - (t880 + t996) * m(3) + t862) * t665 - t728 * t659 + m(3) * t656 * t936 + 0.4e1 * (-t695 - t909 / 0.2e1) * t1002 - (-0.2e1 * (-t689 * t970 + t701 * t731) * t948 - t701 * (t689 * t918 + t701 * t751 + t725 * t915)) * t833 - ((t881 + t906) * t698 + (m(3) * t689 - t701 * t769) * t725 * t919) * t827 - 0.2e1 * (-t695 - t909) * t725;
t879 = (t696 * t817 + (-t978 * t817 + (-t829 * t836 - pkin(1)) * t823) * t720 + t657 * t836 - t663 * t938 - t654) * m(3);
t878 = (t697 * t818 + (-t977 * t818 + (-t831 * t836 - pkin(1)) * t825) * t721 + t658 * t836 - t664 * t937 - t655) * m(3);
t877 = (t698 * t819 + (-t976 * t819 + (-t833 * t836 - pkin(1)) * t827) * t722 + t659 * t836 - t665 * t936 - t656) * m(3);
t1 = [(t924 * t741 + (t877 * t712 - t718 * t999) * t849) * t757 + (t925 * t739 + (-t716 * t1000 + t878 * t710) * t847) * t756 + (t926 * t737 + (-t714 * t1001 + t879 * t708) * t845) * t755; (t924 * t740 + (t877 * t713 + t719 * t999) * t849) * t757 + (t925 * t738 + (t717 * t1000 + t878 * t711) * t847) * t756 + (t926 * t736 + (t715 * t1001 + t879 * t709) * t845) * t755; (t877 * t760 + t921 * t827) * t849 + (t878 * t759 + t922 * t825) * t847 + (t879 * t758 + t923 * t823) * t845;];
taucX  = t1;
