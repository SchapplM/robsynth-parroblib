% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V2G2A0
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:20:33
% EndTime: 2020-08-06 19:20:41
% DurationCPUTime: 8.48s
% Computational Cost: add. (80298->447), mult. (94860->715), div. (11403->6), fcn. (72603->18), ass. (0->295)
t833 = sin(qJ(2,3));
t1028 = 0.2e1 * t833;
t835 = sin(qJ(2,2));
t1027 = 0.2e1 * t835;
t837 = sin(qJ(2,1));
t1026 = 0.2e1 * t837;
t834 = sin(qJ(1,3));
t840 = cos(qJ(1,3));
t852 = pkin(5) - pkin(6);
t879 = pkin(1) * t834 - t852 * t840;
t962 = t834 * qJ(3,3);
t755 = t879 * t833 + t962;
t1003 = pkin(1) * qJ(3,3);
t853 = pkin(2) + pkin(3);
t967 = (qJ(3,3) + t853) * (-qJ(3,3) + t853);
t900 = t833 * t967;
t766 = -t900 + t1003;
t830 = legFrame(3,2);
t799 = sin(t830);
t802 = cos(t830);
t839 = cos(qJ(2,3));
t821 = t839 ^ 2;
t899 = t834 * t967;
t1011 = -0.2e1 * t853;
t915 = qJ(3,3) * t1011;
t788 = pkin(1) * t833 + qJ(3,3);
t948 = t853 * t788;
t973 = t799 * qJ(3,3);
t892 = t833 * t962;
t979 = (t879 + 0.2e1 * t892) * t853;
t719 = (-t799 * t899 + t802 * t915) * t821 + (-t802 * t766 - t799 * t979) * t839 - t755 * t973 + t802 * t948;
t970 = t802 * qJ(3,3);
t720 = (t799 * t915 + t802 * t899) * t821 + (-t799 * t766 + t802 * t979) * t839 + t755 * t970 + t799 * t948;
t791 = t834 * t852;
t796 = t833 * qJ(3,3);
t901 = t821 * t967;
t942 = t853 * t839;
t743 = t840 * t901 + ((0.2e1 * t796 + pkin(1)) * t840 + t791) * t942 + qJ(3,3) * (t788 * t840 + t833 * t791);
t849 = xDP(3);
t850 = xDP(2);
t851 = xDP(1);
t785 = t796 + pkin(1);
t1021 = t785 + t942;
t772 = 0.1e1 / t1021;
t855 = 0.1e1 / qJ(3,3);
t976 = t772 * t855;
t698 = (t719 * t850 + t720 * t851 + t743 * t849) * t976;
t759 = t879 + t892;
t945 = t853 * t833;
t961 = t834 * t853;
t728 = (-t799 * t961 - t970) * t821 + (-t759 * t799 + t802 * t945) * t839 + t802 * t788;
t729 = (t802 * t961 - t973) * t821 + (t759 * t802 + t799 * t945) * t839 + t799 * t788;
t982 = (t1021 * t840 + t791) * t839;
t710 = (t728 * t850 + t729 * t851 + t849 * t982) * t976;
t951 = t853 * t710;
t692 = t698 - t951;
t836 = sin(qJ(1,2));
t842 = cos(qJ(1,2));
t878 = pkin(1) * t836 - t852 * t842;
t958 = t836 * qJ(3,2);
t756 = t878 * t835 + t958;
t1004 = pkin(1) * qJ(3,2);
t966 = (qJ(3,2) + t853) * (-qJ(3,2) + t853);
t897 = t835 * t966;
t767 = -t897 + t1004;
t831 = legFrame(2,2);
t800 = sin(t831);
t803 = cos(t831);
t841 = cos(qJ(2,2));
t822 = t841 ^ 2;
t896 = t836 * t966;
t916 = qJ(3,2) * t1011;
t789 = pkin(1) * t835 + qJ(3,2);
t947 = t853 * t789;
t972 = t800 * qJ(3,2);
t891 = t835 * t958;
t978 = (t878 + 0.2e1 * t891) * t853;
t721 = (-t800 * t896 + t803 * t916) * t822 + (-t803 * t767 - t800 * t978) * t841 - t756 * t972 + t803 * t947;
t969 = t803 * qJ(3,2);
t722 = (t800 * t916 + t803 * t896) * t822 + (-t800 * t767 + t803 * t978) * t841 + t756 * t969 + t800 * t947;
t792 = t836 * t852;
t797 = t835 * qJ(3,2);
t898 = t822 * t966;
t941 = t853 * t841;
t744 = t842 * t898 + ((0.2e1 * t797 + pkin(1)) * t842 + t792) * t941 + qJ(3,2) * (t789 * t842 + t835 * t792);
t786 = t797 + pkin(1);
t1020 = t786 + t941;
t773 = 0.1e1 / t1020;
t857 = 0.1e1 / qJ(3,2);
t975 = t773 * t857;
t699 = (t721 * t850 + t722 * t851 + t744 * t849) * t975;
t761 = t878 + t891;
t944 = t853 * t835;
t957 = t836 * t853;
t730 = (-t800 * t957 - t969) * t822 + (-t761 * t800 + t803 * t944) * t841 + t803 * t789;
t731 = (t803 * t957 - t972) * t822 + (t761 * t803 + t800 * t944) * t841 + t800 * t789;
t981 = (t1020 * t842 + t792) * t841;
t711 = (t730 * t850 + t731 * t851 + t849 * t981) * t975;
t950 = t853 * t711;
t693 = t699 - t950;
t838 = sin(qJ(1,1));
t844 = cos(qJ(1,1));
t877 = pkin(1) * t838 - t852 * t844;
t954 = t838 * qJ(3,1);
t757 = t877 * t837 + t954;
t1005 = pkin(1) * qJ(3,1);
t965 = (qJ(3,1) + t853) * (-qJ(3,1) + t853);
t894 = t837 * t965;
t768 = -t894 + t1005;
t832 = legFrame(1,2);
t801 = sin(t832);
t804 = cos(t832);
t843 = cos(qJ(2,1));
t823 = t843 ^ 2;
t893 = t838 * t965;
t917 = qJ(3,1) * t1011;
t790 = pkin(1) * t837 + qJ(3,1);
t946 = t853 * t790;
t971 = t801 * qJ(3,1);
t890 = t837 * t954;
t977 = (t877 + 0.2e1 * t890) * t853;
t723 = (-t801 * t893 + t804 * t917) * t823 + (-t804 * t768 - t801 * t977) * t843 - t757 * t971 + t804 * t946;
t968 = t804 * qJ(3,1);
t724 = (t801 * t917 + t804 * t893) * t823 + (-t801 * t768 + t804 * t977) * t843 + t757 * t968 + t801 * t946;
t793 = t838 * t852;
t798 = t837 * qJ(3,1);
t895 = t823 * t965;
t940 = t853 * t843;
t745 = t844 * t895 + ((0.2e1 * t798 + pkin(1)) * t844 + t793) * t940 + qJ(3,1) * (t790 * t844 + t837 * t793);
t787 = t798 + pkin(1);
t1019 = t787 + t940;
t774 = 0.1e1 / t1019;
t859 = 0.1e1 / qJ(3,1);
t974 = t774 * t859;
t700 = (t723 * t850 + t724 * t851 + t745 * t849) * t974;
t763 = t877 + t890;
t943 = t853 * t837;
t953 = t838 * t853;
t732 = (-t801 * t953 - t968) * t823 + (-t763 * t801 + t804 * t943) * t843 + t804 * t790;
t733 = (t804 * t953 - t971) * t823 + (t763 * t804 + t801 * t943) * t843 + t801 * t790;
t980 = (t1019 * t844 + t793) * t843;
t712 = (t732 * t850 + t733 * t851 + t849 * t980) * t974;
t949 = t853 * t712;
t694 = t700 - t949;
t1025 = 0.2e1 * t853;
t1015 = 2 * rSges(3,3);
t854 = qJ(3,3) ^ 2;
t1018 = qJ(3,3) * t1015 + t854;
t1009 = t821 - 0.1e1;
t1014 = -0.2e1 * t821;
t939 = t854 * t710;
t683 = t853 * t692 - t939;
t689 = t852 * t692;
t740 = (-t834 * t849 + (-t799 * t850 + t802 * t851) * t840) * t772;
t867 = pkin(1) ^ 2;
t875 = (pkin(6) ^ 2) + t867 + ((-2 * pkin(6) + pkin(5)) * pkin(5));
t775 = t854 + t875;
t826 = t853 ^ 2;
t909 = t740 * t1003;
t925 = 0.2e1 * pkin(1);
t928 = 0.2e1 * t909 + t689;
t952 = t852 * t853;
t963 = t833 * t852;
t989 = t740 * t852;
t990 = t740 * t833;
t991 = t740 * t821;
t904 = t740 * t963;
t994 = (t904 - t951) * t839;
t671 = ((-(t826 - 0.3e1 * t854) * t942 * t991 + (t852 * (-t1025 * t710 + t698) * qJ(3,3) + (-0.3e1 * (-t854 / 0.3e1 + t826) * t796 + (t854 - t826) * t925) * t740) * t821 + (-t939 * t963 + ((-0.4e1 * t909 - t689) * t833 - t740 * (0.3e1 * t854 + t875)) * t853) * t839 - (t775 * t990 + t928) * qJ(3,3)) * t740 + ((t683 * t853 + t900 * t989) * t839 + t683 * pkin(1) + (t683 * t833 + (t1014 + 0.1e1) * t740 * t952) * qJ(3,3)) * t710 + ((pkin(1) * t710 - t994) * t853 + (t1009 * t989 + t710 * t945) * qJ(3,3)) * t698) * t976;
t889 = t839 * qJ(3,3) * t710;
t674 = (-(t852 * t889 + t928 * t833 + (t1025 * t785 * t839 + t775 + t901) * t740) * t839 * t740 + (-qJ(3,3) * t821 * t989 + (t692 + t904) * t942 + t692 * t785) * t710 + (t710 * t785 - t994) * t698) * t976;
t680 = (t692 * t1028 + 0.2e1 * t889 + t989) * t740 * t772;
t827 = -qJ(3,3) - rSges(3,3);
t848 = (m(2) * rSges(2,2));
t887 = -rSges(2,1) * t848 + Icges(2,4) - Icges(3,5);
t846 = (pkin(2) + rSges(3,1));
t995 = t846 * m(3);
t769 = -t827 * t995 + t887;
t704 = t769 * t710;
t707 = t710 ^ 2;
t1010 = (m(2) * rSges(2,1));
t847 = rSges(2,3) + pkin(5);
t845 = pkin(5) + rSges(3,2);
t886 = t845 * t995 - Icges(3,4) - Icges(2,5);
t764 = (t847 * t1010) + t886;
t935 = (Icges(2,6) - Icges(3,6));
t872 = -t847 * t848 + t935;
t998 = t827 * m(3);
t908 = t845 * t998;
t746 = -(t872 - t908) * t839 + t833 * t764;
t861 = rSges(2,2) ^ 2;
t863 = rSges(2,1) ^ 2;
t936 = -Icges(2,1) - Icges(3,1);
t873 = Icges(2,2) + Icges(3,3) + (-t861 + t863) * m(2) + t936;
t919 = rSges(3,3) + t846;
t920 = rSges(3,3) - t846;
t752 = -(qJ(3,3) + t919) * (qJ(3,3) + t920) * m(3) + t873;
t1006 = m(2) * t847;
t765 = (rSges(2,1) * t1006) + t886;
t782 = t995 + t1010;
t999 = t782 * pkin(1);
t778 = 0.2e1 * t999;
t871 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - ((t847 ^ 2) + t861 + t867) * m(2) - Icges(1,3) + t936;
t860 = rSges(3,3) ^ 2;
t884 = (t845 ^ 2) + t860 + t867;
t885 = rSges(2,2) * t1006 - t935;
t914 = t698 * t998;
t1002 = (t848 + t998) * pkin(1);
t923 = 0.2e1 * t1002;
t924 = -0.2e1 * m(3) * t845;
t964 = t833 * t845;
t934 = (-t752 * t821 - (t769 * t1028 + t778) * t839 + t833 * t923 - (t884 + t1018) * m(3) + t871) * t680 + t746 * t674 - m(3) * t671 * t964 - 0.4e1 * (-t704 - t914 / 0.2e1) * t991 + (-0.2e1 * (-t698 * t995 + t752 * t710) * t990 - t710 * (t698 * t924 + t765 * t710 + t740 * t923)) * t839 + ((t885 + t908) * t707 + (m(3) * t698 - t710 * t782) * t740 * t925) * t833 + 0.2e1 * t740 * (-t704 - t914);
t1024 = t840 * t934;
t856 = qJ(3,2) ^ 2;
t1017 = qJ(3,2) * t1015 + t856;
t1008 = t822 - 0.1e1;
t1013 = -0.2e1 * t822;
t938 = t856 * t711;
t684 = t853 * t693 - t938;
t690 = t852 * t693;
t741 = (-t836 * t849 + (-t800 * t850 + t803 * t851) * t842) * t773;
t776 = t856 + t875;
t910 = t741 * t1004;
t927 = 0.2e1 * t910 + t690;
t959 = t835 * t852;
t986 = t741 * t852;
t987 = t741 * t835;
t988 = t741 * t822;
t903 = t741 * t959;
t993 = (t903 - t950) * t841;
t672 = ((-(t826 - 0.3e1 * t856) * t941 * t988 + (t852 * (-t1025 * t711 + t699) * qJ(3,2) + (-0.3e1 * (-t856 / 0.3e1 + t826) * t797 + (t856 - t826) * t925) * t741) * t822 + (-t938 * t959 + ((-0.4e1 * t910 - t690) * t835 - t741 * (0.3e1 * t856 + t875)) * t853) * t841 - (t776 * t987 + t927) * qJ(3,2)) * t741 + ((t684 * t853 + t897 * t986) * t841 + t684 * pkin(1) + (t684 * t835 + (t1013 + 0.1e1) * t741 * t952) * qJ(3,2)) * t711 + ((pkin(1) * t711 - t993) * t853 + (t1008 * t986 + t711 * t944) * qJ(3,2)) * t699) * t975;
t888 = t841 * qJ(3,2) * t711;
t675 = (-(t852 * t888 + t927 * t835 + (t1025 * t786 * t841 + t776 + t898) * t741) * t841 * t741 + (-qJ(3,2) * t822 * t986 + (t693 + t903) * t941 + t693 * t786) * t711 + (t711 * t786 - t993) * t699) * t975;
t681 = (t693 * t1027 + 0.2e1 * t888 + t986) * t741 * t773;
t828 = -qJ(3,2) - rSges(3,3);
t770 = -t828 * t995 + t887;
t705 = t770 * t711;
t708 = t711 ^ 2;
t997 = t828 * m(3);
t907 = t845 * t997;
t747 = -(t872 - t907) * t841 + t835 * t764;
t753 = -(qJ(3,2) + t919) * (qJ(3,2) + t920) * m(3) + t873;
t913 = t699 * t997;
t1001 = (t848 + t997) * pkin(1);
t922 = 0.2e1 * t1001;
t960 = t835 * t845;
t933 = (-t753 * t822 - (t770 * t1027 + t778) * t841 + t835 * t922 - (t884 + t1017) * m(3) + t871) * t681 + t747 * t675 - m(3) * t672 * t960 - 0.4e1 * (-t705 - t913 / 0.2e1) * t988 + (-0.2e1 * (-t699 * t995 + t753 * t711) * t987 - t711 * (t699 * t924 + t765 * t711 + t741 * t922)) * t841 + ((t885 + t907) * t708 + (m(3) * t699 - t711 * t782) * t741 * t925) * t835 + 0.2e1 * t741 * (-t705 - t913);
t1023 = t842 * t933;
t858 = qJ(3,1) ^ 2;
t1016 = qJ(3,1) * t1015 + t858;
t1007 = t823 - 0.1e1;
t1012 = -0.2e1 * t823;
t937 = t858 * t712;
t685 = t853 * t694 - t937;
t691 = t852 * t694;
t742 = (-t838 * t849 + (-t801 * t850 + t804 * t851) * t844) * t774;
t777 = t858 + t875;
t911 = t742 * t1005;
t926 = 0.2e1 * t911 + t691;
t955 = t837 * t852;
t983 = t742 * t852;
t984 = t742 * t837;
t985 = t742 * t823;
t902 = t742 * t955;
t992 = (t902 - t949) * t843;
t673 = ((-(t826 - 0.3e1 * t858) * t940 * t985 + (t852 * (-t1025 * t712 + t700) * qJ(3,1) + (-0.3e1 * (-t858 / 0.3e1 + t826) * t798 + (t858 - t826) * t925) * t742) * t823 + (-t937 * t955 + ((-0.4e1 * t911 - t691) * t837 - t742 * (0.3e1 * t858 + t875)) * t853) * t843 - (t777 * t984 + t926) * qJ(3,1)) * t742 + ((t685 * t853 + t894 * t983) * t843 + t685 * pkin(1) + (t685 * t837 + (t1012 + 0.1e1) * t742 * t952) * qJ(3,1)) * t712 + ((pkin(1) * t712 - t992) * t853 + (t1007 * t983 + t712 * t943) * qJ(3,1)) * t700) * t974;
t905 = t712 * t843 * qJ(3,1);
t676 = (-(t852 * t905 + t926 * t837 + (t1025 * t787 * t843 + t777 + t895) * t742) * t843 * t742 + (-qJ(3,1) * t823 * t983 + (t694 + t902) * t940 + t694 * t787) * t712 + (t712 * t787 - t992) * t700) * t974;
t682 = (t694 * t1026 + 0.2e1 * t905 + t983) * t742 * t774;
t829 = -qJ(3,1) - rSges(3,3);
t771 = -t829 * t995 + t887;
t706 = t771 * t712;
t709 = t712 ^ 2;
t996 = t829 * m(3);
t906 = t845 * t996;
t748 = -(t872 - t906) * t843 + t837 * t764;
t754 = -(qJ(3,1) + t919) * (qJ(3,1) + t920) * m(3) + t873;
t912 = t700 * t996;
t1000 = (t848 + t996) * pkin(1);
t921 = 0.2e1 * t1000;
t956 = t837 * t845;
t932 = (-t754 * t823 - (t771 * t1026 + t778) * t843 + t837 * t921 - (t884 + t1016) * m(3) + t871) * t682 + t748 * t676 - m(3) * t673 * t956 - 0.4e1 * (-t706 - t912 / 0.2e1) * t985 + (-0.2e1 * (-t700 * t995 + t754 * t712) * t984 - t712 * (t700 * t924 + t765 * t712 + t742 * t921)) * t843 + ((t885 + t906) * t709 + (m(3) * t700 - t712 * t782) * t742 * t925) * t837 + 0.2e1 * t742 * (-t706 - t912);
t1022 = t844 * t932;
t737 = t740 ^ 2;
t874 = -(t861 + t863) * m(2) - Icges(3,2) - Icges(2,3);
t876 = pkin(2) ^ 2 + t860 + (2 * pkin(2) + rSges(3,1)) * rSges(3,1);
t931 = t746 * t680 + (-(t876 + t1018) * m(3) + t874) * t674 + t671 * t995 - 0.2e1 * t710 * t914 + (t769 * t1014 + (t752 * t833 + t1002) * t839 + t833 * t999 + t769) * t737;
t738 = t741 ^ 2;
t930 = t747 * t681 + (-(t876 + t1017) * m(3) + t874) * t675 + t672 * t995 - 0.2e1 * t711 * t913 + (t770 * t1013 + (t753 * t835 + t1001) * t841 + t835 * t999 + t770) * t738;
t739 = t742 ^ 2;
t929 = t748 * t682 + (-(t876 + t1016) * m(3) + t874) * t676 + t673 * t995 - 0.2e1 * t712 * t912 + (t771 * t1012 + (t754 * t837 + t1000) * t843 + t837 * t999 + t771) * t739;
t883 = (t707 * t827 + (-t1009 * t827 + (-t839 * t846 - pkin(1)) * t833) * t737 + t674 * t846 - t680 * t964 - t671) * m(3);
t882 = (t708 * t828 + (-t1008 * t828 + (-t841 * t846 - pkin(1)) * t835) * t738 + t675 * t846 - t681 * t960 - t672) * m(3);
t881 = (t709 * t829 + (-t1007 * t829 + (-t843 * t846 - pkin(1)) * t837) * t739 + t676 * t846 - t682 * t956 - t673) * m(3);
t1 = [(t804 * t1022 + (t881 * t724 + t929 * t733) * t859) * t774 + (t803 * t1023 + (t882 * t722 + t930 * t731) * t857) * t773 + (t802 * t1024 + (t883 * t720 + t931 * t729) * t855) * t772; (-t801 * t1022 + (t881 * t723 + t929 * t732) * t859) * t774 + (-t800 * t1023 + (t882 * t721 + t930 * t730) * t857) * t773 + (-t799 * t1024 + (t883 * t719 + t931 * t728) * t855) * t772; (-t932 * t838 + (t881 * t745 + t929 * t980) * t859) * t774 + (-t933 * t836 + (t882 * t744 + t930 * t981) * t857) * t773 + (-t934 * t834 + (t883 * t743 + t931 * t982) * t855) * t772;];
taucX  = t1;
