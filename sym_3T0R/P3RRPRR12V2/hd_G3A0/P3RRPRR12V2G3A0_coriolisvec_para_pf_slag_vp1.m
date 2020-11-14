% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V2G3A0
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
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:27:00
% EndTime: 2020-08-06 19:27:08
% DurationCPUTime: 8.69s
% Computational Cost: add. (80298->447), mult. (94860->715), div. (11403->6), fcn. (72603->18), ass. (0->298)
t833 = sin(qJ(2,3));
t1031 = 0.2e1 * t833;
t835 = sin(qJ(2,2));
t1030 = 0.2e1 * t835;
t837 = sin(qJ(2,1));
t1029 = 0.2e1 * t837;
t840 = cos(qJ(1,3));
t834 = sin(qJ(1,3));
t852 = pkin(5) - pkin(6);
t925 = pkin(1) * t840 + t834 * t852;
t755 = qJ(3,3) * t840 + t925 * t833;
t796 = t833 * qJ(3,3);
t900 = t840 * t796;
t758 = 0.2e1 * t900 + t925;
t1001 = pkin(1) * qJ(3,3);
t853 = pkin(2) + pkin(3);
t958 = (qJ(3,3) + t853) * (-qJ(3,3) + t853);
t894 = t833 * t958;
t766 = -t894 + t1001;
t997 = pkin(1) * t833;
t788 = qJ(3,3) + t997;
t830 = legFrame(3,2);
t799 = sin(t830);
t802 = cos(t830);
t839 = cos(qJ(2,3));
t821 = t839 ^ 2;
t893 = t840 * t958;
t1014 = -0.2e1 * t853;
t912 = qJ(3,3) * t1014;
t961 = t802 * t853;
t964 = t799 * t853;
t990 = qJ(3,3) * t799;
t719 = (-t799 * t893 + t802 * t912) * t821 + (-t758 * t964 - t766 * t802) * t839 - t755 * t990 + t788 * t961;
t989 = qJ(3,3) * t802;
t720 = (t799 * t912 + t802 * t893) * t821 + (t758 * t961 - t766 * t799) * t839 + t755 * t989 + t788 * t964;
t895 = t821 * t958;
t940 = t852 * t840;
t946 = t839 * t853;
t743 = -t834 * t895 - ((0.2e1 * t796 + pkin(1)) * t834 - t940) * t946 - qJ(3,3) * (t788 * t834 - t833 * t940);
t849 = xDP(3);
t850 = xDP(2);
t851 = xDP(1);
t785 = t796 + pkin(1);
t1024 = t785 + t946;
t772 = 0.1e1 / t1024;
t855 = 0.1e1 / qJ(3,3);
t967 = t772 * t855;
t698 = (t719 * t850 + t720 * t851 + t743 * t849) * t967;
t759 = t900 + t925;
t945 = t840 * t853;
t953 = t833 * t853;
t728 = (-t799 * t945 - t989) * t821 + (-t759 * t799 + t802 * t953) * t839 + t802 * t788;
t729 = (t802 * t945 - t990) * t821 + (t759 * t802 + t799 * t953) * t839 + t799 * t788;
t970 = (t1024 * t834 - t940) * t839;
t710 = (t728 * t850 + t729 * t851 - t849 * t970) * t967;
t985 = t710 * t853;
t692 = t698 - t985;
t842 = cos(qJ(1,2));
t836 = sin(qJ(1,2));
t924 = pkin(1) * t842 + t836 * t852;
t756 = qJ(3,2) * t842 + t924 * t835;
t797 = t835 * qJ(3,2);
t901 = t842 * t797;
t760 = 0.2e1 * t901 + t924;
t1002 = pkin(1) * qJ(3,2);
t957 = (qJ(3,2) + t853) * (-qJ(3,2) + t853);
t891 = t835 * t957;
t767 = -t891 + t1002;
t996 = pkin(1) * t835;
t789 = qJ(3,2) + t996;
t831 = legFrame(2,2);
t800 = sin(t831);
t803 = cos(t831);
t841 = cos(qJ(2,2));
t822 = t841 ^ 2;
t890 = t842 * t957;
t913 = qJ(3,2) * t1014;
t960 = t803 * t853;
t963 = t800 * t853;
t992 = qJ(3,2) * t800;
t721 = (-t800 * t890 + t803 * t913) * t822 + (-t760 * t963 - t767 * t803) * t841 - t756 * t992 + t789 * t960;
t991 = qJ(3,2) * t803;
t722 = (t800 * t913 + t803 * t890) * t822 + (t760 * t960 - t767 * t800) * t841 + t756 * t991 + t789 * t963;
t892 = t822 * t957;
t939 = t852 * t842;
t944 = t841 * t853;
t744 = -t836 * t892 - ((0.2e1 * t797 + pkin(1)) * t836 - t939) * t944 - qJ(3,2) * (t789 * t836 - t835 * t939);
t786 = t797 + pkin(1);
t1023 = t786 + t944;
t773 = 0.1e1 / t1023;
t857 = 0.1e1 / qJ(3,2);
t966 = t773 * t857;
t699 = (t721 * t850 + t722 * t851 + t744 * t849) * t966;
t761 = t901 + t924;
t943 = t842 * t853;
t950 = t835 * t853;
t730 = (-t800 * t943 - t991) * t822 + (-t761 * t800 + t803 * t950) * t841 + t803 * t789;
t731 = (t803 * t943 - t992) * t822 + (t761 * t803 + t800 * t950) * t841 + t800 * t789;
t969 = (t1023 * t836 - t939) * t841;
t711 = (t730 * t850 + t731 * t851 - t849 * t969) * t966;
t983 = t711 * t853;
t693 = t699 - t983;
t844 = cos(qJ(1,1));
t838 = sin(qJ(1,1));
t923 = pkin(1) * t844 + t838 * t852;
t757 = qJ(3,1) * t844 + t923 * t837;
t798 = t837 * qJ(3,1);
t902 = t844 * t798;
t762 = 0.2e1 * t902 + t923;
t1003 = pkin(1) * qJ(3,1);
t956 = (qJ(3,1) + t853) * (-qJ(3,1) + t853);
t888 = t837 * t956;
t768 = -t888 + t1003;
t995 = pkin(1) * t837;
t790 = qJ(3,1) + t995;
t832 = legFrame(1,2);
t801 = sin(t832);
t804 = cos(t832);
t843 = cos(qJ(2,1));
t823 = t843 ^ 2;
t887 = t844 * t956;
t914 = qJ(3,1) * t1014;
t959 = t804 * t853;
t962 = t801 * t853;
t994 = qJ(3,1) * t801;
t723 = (-t801 * t887 + t804 * t914) * t823 + (-t762 * t962 - t768 * t804) * t843 - t757 * t994 + t790 * t959;
t993 = qJ(3,1) * t804;
t724 = (t801 * t914 + t804 * t887) * t823 + (t762 * t959 - t768 * t801) * t843 + t757 * t993 + t790 * t962;
t889 = t823 * t956;
t938 = t852 * t844;
t942 = t843 * t853;
t745 = -t838 * t889 - ((0.2e1 * t798 + pkin(1)) * t838 - t938) * t942 - qJ(3,1) * (t790 * t838 - t837 * t938);
t787 = t798 + pkin(1);
t1022 = t787 + t942;
t774 = 0.1e1 / t1022;
t859 = 0.1e1 / qJ(3,1);
t965 = t774 * t859;
t700 = (t723 * t850 + t724 * t851 + t745 * t849) * t965;
t763 = t902 + t923;
t941 = t844 * t853;
t947 = t837 * t853;
t732 = (-t801 * t941 - t993) * t823 + (-t763 * t801 + t804 * t947) * t843 + t804 * t790;
t733 = (t804 * t941 - t994) * t823 + (t763 * t804 + t801 * t947) * t843 + t801 * t790;
t968 = (t1022 * t838 - t938) * t843;
t712 = (t732 * t850 + t733 * t851 - t849 * t968) * t965;
t981 = t712 * t853;
t694 = t700 - t981;
t922 = 0.2e1 * pkin(1);
t1028 = 0.2e1 * t853;
t846 = (pkin(2) + rSges(3,1));
t1004 = m(3) * t846;
t1018 = 2 * rSges(3,3);
t854 = qJ(3,3) ^ 2;
t1021 = qJ(3,3) * t1018 + t854;
t1012 = t821 - 0.1e1;
t1017 = -0.2e1 * t821;
t984 = t710 * t854;
t683 = t692 * t853 - t984;
t689 = t852 * t692;
t740 = (-t840 * t849 + (t799 * t850 - t802 * t851) * t834) * t772;
t867 = pkin(1) ^ 2;
t875 = (pkin(6) ^ 2) + t867 + ((-2 * pkin(6) + pkin(5)) * pkin(5));
t775 = t854 + t875;
t826 = t853 ^ 2;
t903 = t740 * t1001;
t928 = 0.2e1 * t903 + t689;
t937 = t852 * t853;
t954 = t833 * t852;
t977 = t740 * t852;
t978 = t740 * t833;
t979 = t740 * t821;
t898 = t740 * t954;
t988 = (t898 - t985) * t839;
t671 = ((-(t826 - 0.3e1 * t854) * t946 * t979 + (t852 * (-t1028 * t710 + t698) * qJ(3,3) + (-0.3e1 * (-t854 / 0.3e1 + t826) * t796 + (t854 - t826) * t922) * t740) * t821 + (-t954 * t984 + ((-0.4e1 * t903 - t689) * t833 - t740 * (0.3e1 * t854 + t875)) * t853) * t839 - qJ(3,3) * (t775 * t978 + t928)) * t740 + ((t683 * t853 + t894 * t977) * t839 + t683 * pkin(1) + (t683 * t833 + (t1017 + 0.1e1) * t740 * t937) * qJ(3,3)) * t710 + ((pkin(1) * t710 - t988) * t853 + (t1012 * t977 + t710 * t953) * qJ(3,3)) * t698) * t967;
t886 = t839 * qJ(3,3) * t710;
t674 = (-(t852 * t886 + t928 * t833 + (t1028 * t785 * t839 + t775 + t895) * t740) * t839 * t740 + (-qJ(3,3) * t821 * t977 + (t692 + t898) * t946 + t692 * t785) * t710 + (t710 * t785 - t988) * t698) * t967;
t680 = (t692 * t1031 + 0.2e1 * t886 + t977) * t740 * t772;
t827 = -qJ(3,3) - rSges(3,3);
t848 = (m(2) * rSges(2,2));
t884 = -rSges(2,1) * t848 + Icges(2,4) - Icges(3,5);
t769 = -t827 * t1004 + t884;
t704 = t769 * t710;
t707 = t710 ^ 2;
t1013 = (m(2) * rSges(2,1));
t847 = rSges(2,3) + pkin(5);
t845 = pkin(5) + rSges(3,2);
t883 = t845 * t1004 - Icges(3,4) - Icges(2,5);
t764 = (t847 * t1013) + t883;
t935 = (Icges(2,6) - Icges(3,6));
t872 = -t847 * t848 + t935;
t1005 = m(3) * t845;
t908 = t827 * t1005;
t746 = -(t872 - t908) * t839 + t833 * t764;
t861 = rSges(2,2) ^ 2;
t863 = rSges(2,1) ^ 2;
t936 = -Icges(2,1) - Icges(3,1);
t873 = Icges(2,2) + Icges(3,3) + (-t861 + t863) * m(2) + t936;
t916 = rSges(3,3) + t846;
t917 = rSges(3,3) - t846;
t752 = -(qJ(3,3) + t916) * (qJ(3,3) + t917) * m(3) + t873;
t1009 = m(2) * t847;
t765 = (rSges(2,1) * t1009) + t883;
t782 = t1004 + t1013;
t778 = t782 * t922;
t871 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - ((t847 ^ 2) + t861 + t867) * m(2) - Icges(1,3) + t936;
t860 = rSges(3,3) ^ 2;
t881 = (t845 ^ 2) + t860 + t867;
t882 = rSges(2,2) * t1009 - t935;
t1008 = m(3) * t827;
t911 = t698 * t1008;
t1000 = pkin(1) * (t848 + t1008);
t920 = 0.2e1 * t1000;
t921 = -2 * t1005;
t955 = t833 * t845;
t934 = (-t752 * t821 - (t769 * t1031 + t778) * t839 + t833 * t920 - (t881 + t1021) * m(3) + t871) * t680 + t746 * t674 - m(3) * t671 * t955 - 0.4e1 * (-t704 - t911 / 0.2e1) * t979 + (-0.2e1 * (-t698 * t1004 + t710 * t752) * t978 - t710 * (t698 * t921 + t710 * t765 + t740 * t920)) * t839 + ((t882 + t908) * t707 + (m(3) * t698 - t710 * t782) * t740 * t922) * t833 + 0.2e1 * t740 * (-t704 - t911);
t1027 = t834 * t934;
t856 = qJ(3,2) ^ 2;
t1020 = qJ(3,2) * t1018 + t856;
t1011 = t822 - 0.1e1;
t1016 = -0.2e1 * t822;
t982 = t711 * t856;
t684 = t853 * t693 - t982;
t690 = t852 * t693;
t741 = (-t842 * t849 + (t800 * t850 - t803 * t851) * t836) * t773;
t776 = t856 + t875;
t904 = t741 * t1002;
t927 = 0.2e1 * t904 + t690;
t951 = t835 * t852;
t974 = t741 * t852;
t975 = t741 * t835;
t976 = t741 * t822;
t897 = t741 * t951;
t987 = (t897 - t983) * t841;
t672 = ((-(t826 - 0.3e1 * t856) * t944 * t976 + (t852 * (-t1028 * t711 + t699) * qJ(3,2) + (-0.3e1 * (-t856 / 0.3e1 + t826) * t797 + (t856 - t826) * t922) * t741) * t822 + (-t951 * t982 + ((-0.4e1 * t904 - t690) * t835 - t741 * (0.3e1 * t856 + t875)) * t853) * t841 - qJ(3,2) * (t776 * t975 + t927)) * t741 + ((t684 * t853 + t891 * t974) * t841 + t684 * pkin(1) + (t684 * t835 + (t1016 + 0.1e1) * t741 * t937) * qJ(3,2)) * t711 + ((pkin(1) * t711 - t987) * t853 + (t1011 * t974 + t711 * t950) * qJ(3,2)) * t699) * t966;
t885 = t841 * qJ(3,2) * t711;
t675 = (-(t852 * t885 + t927 * t835 + (t1028 * t786 * t841 + t776 + t892) * t741) * t841 * t741 + (-qJ(3,2) * t822 * t974 + (t693 + t897) * t944 + t693 * t786) * t711 + (t711 * t786 - t987) * t699) * t966;
t681 = (t693 * t1030 + 0.2e1 * t885 + t974) * t741 * t773;
t828 = -qJ(3,2) - rSges(3,3);
t770 = -t828 * t1004 + t884;
t705 = t770 * t711;
t708 = t711 ^ 2;
t907 = t828 * t1005;
t747 = -(t872 - t907) * t841 + t835 * t764;
t753 = -(qJ(3,2) + t916) * (qJ(3,2) + t917) * m(3) + t873;
t1007 = m(3) * t828;
t910 = t699 * t1007;
t999 = pkin(1) * (t848 + t1007);
t919 = 0.2e1 * t999;
t952 = t835 * t845;
t933 = (-t753 * t822 - (t770 * t1030 + t778) * t841 + t835 * t919 - (t881 + t1020) * m(3) + t871) * t681 + t747 * t675 - m(3) * t672 * t952 - 0.4e1 * (-t705 - t910 / 0.2e1) * t976 + (-0.2e1 * (-t699 * t1004 + t711 * t753) * t975 - t711 * (t699 * t921 + t711 * t765 + t741 * t919)) * t841 + ((t882 + t907) * t708 + (m(3) * t699 - t711 * t782) * t741 * t922) * t835 + 0.2e1 * t741 * (-t705 - t910);
t1026 = t836 * t933;
t858 = qJ(3,1) ^ 2;
t1019 = qJ(3,1) * t1018 + t858;
t1010 = t823 - 0.1e1;
t1015 = -0.2e1 * t823;
t980 = t712 * t858;
t685 = t853 * t694 - t980;
t691 = t852 * t694;
t742 = (-t844 * t849 + (t801 * t850 - t804 * t851) * t838) * t774;
t777 = t858 + t875;
t905 = t742 * t1003;
t926 = 0.2e1 * t905 + t691;
t948 = t837 * t852;
t971 = t742 * t852;
t972 = t742 * t837;
t973 = t742 * t823;
t896 = t742 * t948;
t986 = (t896 - t981) * t843;
t673 = ((-(t826 - 0.3e1 * t858) * t942 * t973 + (t852 * (-t1028 * t712 + t700) * qJ(3,1) + (-0.3e1 * (-t858 / 0.3e1 + t826) * t798 + (t858 - t826) * t922) * t742) * t823 + (-t948 * t980 + ((-0.4e1 * t905 - t691) * t837 - t742 * (0.3e1 * t858 + t875)) * t853) * t843 - qJ(3,1) * (t777 * t972 + t926)) * t742 + ((t685 * t853 + t888 * t971) * t843 + t685 * pkin(1) + (t685 * t837 + (t1015 + 0.1e1) * t742 * t937) * qJ(3,1)) * t712 + ((pkin(1) * t712 - t986) * t853 + (t1010 * t971 + t712 * t947) * qJ(3,1)) * t700) * t965;
t899 = t712 * t843 * qJ(3,1);
t676 = (-(t852 * t899 + t926 * t837 + (t1028 * t787 * t843 + t777 + t889) * t742) * t843 * t742 + (-qJ(3,1) * t823 * t971 + (t694 + t896) * t942 + t694 * t787) * t712 + (t712 * t787 - t986) * t700) * t965;
t682 = (t694 * t1029 + 0.2e1 * t899 + t971) * t742 * t774;
t829 = -qJ(3,1) - rSges(3,3);
t771 = -t829 * t1004 + t884;
t706 = t771 * t712;
t709 = t712 ^ 2;
t906 = t829 * t1005;
t748 = -(t872 - t906) * t843 + t837 * t764;
t754 = -(qJ(3,1) + t916) * (qJ(3,1) + t917) * m(3) + t873;
t1006 = m(3) * t829;
t909 = t700 * t1006;
t998 = pkin(1) * (t848 + t1006);
t918 = 0.2e1 * t998;
t949 = t837 * t845;
t932 = (-t754 * t823 - (t771 * t1029 + t778) * t843 + t837 * t918 - (t881 + t1019) * m(3) + t871) * t682 + t748 * t676 - m(3) * t673 * t949 - 0.4e1 * (-t706 - t909 / 0.2e1) * t973 + (-0.2e1 * (-t700 * t1004 + t712 * t754) * t972 - t712 * (t700 * t921 + t712 * t765 + t742 * t918)) * t843 + ((t882 + t906) * t709 + (m(3) * t700 - t712 * t782) * t742 * t922) * t837 + 0.2e1 * t742 * (-t706 - t909);
t1025 = t838 * t932;
t737 = t740 ^ 2;
t874 = -(t861 + t863) * m(2) - Icges(3,2) - Icges(2,3);
t876 = pkin(2) ^ 2 + t860 + (2 * pkin(2) + rSges(3,1)) * rSges(3,1);
t931 = t746 * t680 + (-(t876 + t1021) * m(3) + t874) * t674 + t671 * t1004 - 0.2e1 * t710 * t911 + (t769 * t1017 + (t752 * t833 + t1000) * t839 + t782 * t997 + t769) * t737;
t738 = t741 ^ 2;
t930 = t747 * t681 + (-(t876 + t1020) * m(3) + t874) * t675 + t672 * t1004 - 0.2e1 * t711 * t910 + (t770 * t1016 + (t753 * t835 + t999) * t841 + t782 * t996 + t770) * t738;
t739 = t742 ^ 2;
t929 = t748 * t682 + (-(t876 + t1019) * m(3) + t874) * t676 + t673 * t1004 - 0.2e1 * t712 * t909 + (t771 * t1015 + (t754 * t837 + t998) * t843 + t782 * t995 + t771) * t739;
t880 = (t707 * t827 + (-t1012 * t827 + (-t839 * t846 - pkin(1)) * t833) * t737 + t674 * t846 - t680 * t955 - t671) * m(3);
t879 = (t708 * t828 + (-t1011 * t828 + (-t841 * t846 - pkin(1)) * t835) * t738 + t675 * t846 - t681 * t952 - t672) * m(3);
t878 = (t709 * t829 + (-t1010 * t829 + (-t843 * t846 - pkin(1)) * t837) * t739 + t676 * t846 - t682 * t949 - t673) * m(3);
t1 = [(-t804 * t1025 + (t878 * t724 + t929 * t733) * t859) * t774 + (-t803 * t1026 + (t879 * t722 + t930 * t731) * t857) * t773 + (-t802 * t1027 + (t880 * t720 + t931 * t729) * t855) * t772; (t801 * t1025 + (t878 * t723 + t929 * t732) * t859) * t774 + (t800 * t1026 + (t879 * t721 + t930 * t730) * t857) * t773 + (t799 * t1027 + (t880 * t719 + t931 * t728) * t855) * t772; (-t932 * t844 + (t878 * t745 - t929 * t968) * t859) * t774 + (-t933 * t842 + (t879 * t744 - t930 * t969) * t857) * t773 + (-t934 * t840 + (t880 * t743 - t931 * t970) * t855) * t772;];
taucX  = t1;
