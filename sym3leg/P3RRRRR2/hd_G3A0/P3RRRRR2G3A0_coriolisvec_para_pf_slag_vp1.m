% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRRRR2G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:11:46
% EndTime: 2020-03-09 21:11:49
% DurationCPUTime: 3.52s
% Computational Cost: add. (9738->276), mult. (24873->548), div. (9648->11), fcn. (27072->45), ass. (0->272)
t845 = legFrame(1,2);
t815 = sin(t845);
t818 = cos(t845);
t861 = cos(qJ(3,1));
t852 = sin(qJ(3,1));
t862 = cos(qJ(2,1));
t993 = t862 * pkin(1);
t933 = t852 * t993;
t853 = sin(qJ(2,1));
t854 = sin(qJ(1,1));
t863 = cos(qJ(1,1));
t784 = t854 * t853 - t863 * t862;
t834 = t861 ^ 2;
t936 = pkin(2) * t784 * t834;
t960 = t818 * t852;
t999 = pkin(1) * t863;
t739 = -t815 * t936 + (-pkin(2) * t960 + t815 * t999) * t861 - t818 * t933;
t869 = xDP(2);
t836 = 0.1e1 / t861 ^ 2;
t827 = 0.1e1 / t853;
t878 = 0.1e1 / pkin(2);
t880 = 0.1e1 / pkin(1);
t946 = t878 * t880;
t918 = t827 * t946;
t902 = t836 * t918;
t733 = t739 * t869 * t902;
t963 = t815 * t852;
t742 = t818 * t936 + (-pkin(2) * t963 - t818 * t999) * t861 - t815 * t933;
t870 = xDP(1);
t736 = t742 * t870 * t902;
t769 = pkin(2) * (t863 * t853 + t854 * t862) * t861 + t854 * pkin(1);
t868 = xDP(3);
t835 = 0.1e1 / t861;
t903 = t835 * t918;
t751 = t769 * t868 * t903;
t721 = t736 + t733 + t751;
t975 = t784 * t861;
t759 = t815 * t975 + t960;
t760 = -t818 * t975 + t963;
t957 = t827 * t880;
t919 = t835 * t957;
t947 = t868 * t880;
t966 = sin(qJ(1,1) + qJ(2,1)) * t827;
t727 = -t947 * t966 + (t759 * t869 + t760 * t870) * t919;
t940 = -t721 - t727;
t1018 = t940 ^ 2;
t860 = cos(qJ(1,2));
t1000 = pkin(1) * t860;
t844 = legFrame(2,2);
t814 = sin(t844);
t817 = cos(t844);
t858 = cos(qJ(3,2));
t849 = sin(qJ(3,2));
t859 = cos(qJ(2,2));
t994 = t859 * pkin(1);
t934 = t849 * t994;
t850 = sin(qJ(2,2));
t851 = sin(qJ(1,2));
t783 = t851 * t850 - t860 * t859;
t831 = t858 ^ 2;
t937 = pkin(2) * t783 * t831;
t961 = t817 * t849;
t738 = -t814 * t937 + (-pkin(2) * t961 + t814 * t1000) * t858 - t817 * t934;
t833 = 0.1e1 / t858 ^ 2;
t826 = 0.1e1 / t850;
t920 = t826 * t946;
t904 = t833 * t920;
t732 = t738 * t869 * t904;
t964 = t814 * t849;
t741 = t817 * t937 + (-pkin(2) * t964 - t817 * t1000) * t858 - t814 * t934;
t735 = t741 * t870 * t904;
t768 = pkin(2) * (t860 * t850 + t851 * t859) * t858 + t851 * pkin(1);
t832 = 0.1e1 / t858;
t905 = t832 * t920;
t750 = t768 * t868 * t905;
t720 = t735 + t732 + t750;
t976 = t783 * t858;
t757 = t814 * t976 + t961;
t758 = -t817 * t976 + t964;
t958 = t826 * t880;
t921 = t832 * t958;
t967 = sin(qJ(1,2) + qJ(2,2)) * t826;
t726 = -t947 * t967 + (t757 * t869 + t758 * t870) * t921;
t941 = -t720 - t726;
t1017 = t941 ^ 2;
t857 = cos(qJ(1,3));
t1001 = pkin(1) * t857;
t843 = legFrame(3,2);
t813 = sin(t843);
t816 = cos(t843);
t855 = cos(qJ(3,3));
t846 = sin(qJ(3,3));
t856 = cos(qJ(2,3));
t995 = t856 * pkin(1);
t935 = t846 * t995;
t847 = sin(qJ(2,3));
t848 = sin(qJ(1,3));
t782 = t848 * t847 - t857 * t856;
t828 = t855 ^ 2;
t938 = pkin(2) * t782 * t828;
t962 = t816 * t846;
t737 = -t813 * t938 + (-pkin(2) * t962 + t813 * t1001) * t855 - t816 * t935;
t830 = 0.1e1 / t855 ^ 2;
t825 = 0.1e1 / t847;
t922 = t825 * t946;
t906 = t830 * t922;
t731 = t737 * t869 * t906;
t965 = t813 * t846;
t740 = t816 * t938 + (-pkin(2) * t965 - t816 * t1001) * t855 - t813 * t935;
t734 = t740 * t870 * t906;
t767 = pkin(2) * (t857 * t847 + t848 * t856) * t855 + t848 * pkin(1);
t829 = 0.1e1 / t855;
t907 = t829 * t922;
t749 = t767 * t868 * t907;
t719 = t734 + t731 + t749;
t977 = t782 * t855;
t755 = t813 * t977 + t962;
t756 = -t816 * t977 + t965;
t959 = t825 * t880;
t923 = t829 * t959;
t968 = sin(qJ(1,3) + qJ(2,3)) * t825;
t725 = -t947 * t968 + (t755 * t869 + t756 * t870) * t923;
t942 = -t719 - t725;
t1016 = t942 ^ 2;
t1015 = 0.2e1 * pkin(1);
t875 = rSges(3,2) ^ 2;
t876 = rSges(3,1) ^ 2;
t788 = (-t875 + t876) * m(3) - Icges(3,1) + Icges(3,2);
t1008 = -t788 / 0.2e1;
t797 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t872 = 0.2e1 * qJ(3,3);
t819 = sin(t872);
t822 = cos(t872);
t1014 = t822 * t1008 + t797 * t819;
t873 = 0.2e1 * qJ(3,2);
t820 = sin(t873);
t823 = cos(t873);
t1013 = t823 * t1008 + t797 * t820;
t874 = 0.2e1 * qJ(3,1);
t821 = sin(t874);
t824 = cos(t874);
t1012 = t824 * t1008 + t797 * t821;
t722 = t725 ^ 2;
t723 = t726 ^ 2;
t724 = t727 ^ 2;
t1011 = -0.2e1 * t788;
t1010 = -m(3) / 0.2e1;
t867 = m(2) * rSges(2,1);
t1009 = m(3) * rSges(3,3);
t837 = qJ(3,3) + qJ(2,3);
t1007 = sin(t837) / 0.2e1;
t839 = qJ(2,2) + qJ(3,2);
t1006 = sin(t839) / 0.2e1;
t841 = qJ(3,1) + qJ(2,1);
t1005 = sin(t841) / 0.2e1;
t1004 = cos(t837) / 0.2e1;
t1003 = cos(t839) / 0.2e1;
t1002 = cos(t841) / 0.2e1;
t998 = t722 * pkin(1);
t997 = t723 * pkin(1);
t996 = t724 * pkin(1);
t764 = (t813 * t870 + t816 * t869) * t878 * t829;
t761 = t764 ^ 2;
t710 = t734 / 0.2e1 + t731 / 0.2e1 + t749 / 0.2e1 + t725;
t986 = t719 * t710;
t992 = (t761 / 0.2e1 + t986) * t847;
t765 = (t814 * t870 + t817 * t869) * t878 * t832;
t762 = t765 ^ 2;
t711 = t735 / 0.2e1 + t732 / 0.2e1 + t750 / 0.2e1 + t726;
t985 = t720 * t711;
t991 = (t762 / 0.2e1 + t985) * t850;
t766 = (t815 * t870 + t818 * t869) * t878 * t835;
t763 = t766 ^ 2;
t712 = t736 / 0.2e1 + t733 / 0.2e1 + t751 / 0.2e1 + t727;
t984 = t721 * t712;
t990 = (t763 / 0.2e1 + t984) * t853;
t989 = t942 * t828;
t988 = t941 * t831;
t987 = t940 * t834;
t983 = t761 * t829;
t982 = t762 * t832;
t981 = t763 * t835;
t980 = t764 * t942;
t979 = t765 * t941;
t978 = t766 * t940;
t974 = t788 * t819;
t973 = t788 * t820;
t972 = t788 * t821;
t971 = t797 * t822;
t970 = t797 * t823;
t969 = t797 * t824;
t956 = t846 * t847;
t955 = t849 * t850;
t954 = t852 * t853;
t953 = t855 * t856;
t952 = t856 * t764;
t951 = t858 * t859;
t950 = t859 * t765;
t949 = t861 * t862;
t948 = t862 * t766;
t877 = pkin(2) ^ 2;
t926 = t764 * t956;
t695 = (-t877 * t989 + (pkin(1) * t725 + (0.2e1 * t710 * t953 - t926) * pkin(2)) * pkin(1)) * t725 * t907 + (-pkin(2) * t989 + (-t942 * t953 - t926) * pkin(1)) * t719 * t923 + ((pkin(1) * t942 * t956 + pkin(2) * t764) * t855 + pkin(1) * t952) * t830 * t764 * t959;
t701 = (-t722 * t995 + (-t855 * t1016 - t983) * pkin(2)) * t959;
t901 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t939 = t875 + t876;
t917 = 0.2e1 * rSges(3,3) ^ 2 + t939;
t894 = t917 * t1010 + t901;
t884 = t894 + t1014;
t791 = m(2) * rSges(2,2) - t1009;
t890 = ((-t867 + (-rSges(3,1) * t855 + rSges(3,2) * t846) * m(3)) * t856 + t791 * t847) * pkin(1);
t728 = t890 + t884;
t752 = -t846 * Icges(3,5) - Icges(3,6) * t855 + (rSges(3,1) * t846 + rSges(3,2) * t855) * m(3) * (pkin(1) * t847 + rSges(3,3));
t879 = pkin(1) ^ 2;
t881 = -m(2) * t879 - Icges(1,3) + (0.2e1 * t879 + t917) * t1010 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t901;
t795 = -rSges(3,2) * t1009 + Icges(3,6);
t796 = rSges(3,1) * t1009 - Icges(3,5);
t887 = (-t795 * t846 - t796 * t855) * t761;
t893 = -t791 * t856 - t847 * t867;
t929 = t846 * t983;
t932 = t942 * t952;
t945 = t887 - (t846 * t855 * t1011 + (-0.4e1 * t828 + 0.2e1) * t797) * t980 + (t893 * t986 + ((-rSges(3,1) * t992 + rSges(3,2) * t932) * t855 + (rSges(3,1) * t932 + rSges(3,2) * t992) * t846) * m(3)) * t1015 + (t881 + 0.2e1 * t890 + t1014) * t701 + t728 * t695 - t752 * t929;
t925 = t765 * t955;
t696 = (-t877 * t988 + (pkin(1) * t726 + (0.2e1 * t711 * t951 - t925) * pkin(2)) * pkin(1)) * t726 * t905 + (-pkin(2) * t988 + (-t941 * t951 - t925) * pkin(1)) * t720 * t921 + ((pkin(1) * t941 * t955 + pkin(2) * t765) * t858 + pkin(1) * t950) * t833 * t765 * t958;
t702 = (-t723 * t994 + (-t858 * t1017 - t982) * pkin(2)) * t958;
t883 = t894 + t1013;
t889 = ((-t867 + (-rSges(3,1) * t858 + rSges(3,2) * t849) * m(3)) * t859 + t791 * t850) * pkin(1);
t729 = t889 + t883;
t753 = -t849 * Icges(3,5) - Icges(3,6) * t858 + (rSges(3,1) * t849 + rSges(3,2) * t858) * m(3) * (pkin(1) * t850 + rSges(3,3));
t886 = (-t795 * t849 - t796 * t858) * t762;
t892 = -t791 * t859 - t850 * t867;
t928 = t849 * t982;
t931 = t941 * t950;
t944 = t886 - (t849 * t858 * t1011 + (-0.4e1 * t831 + 0.2e1) * t797) * t979 + (t892 * t985 + ((-rSges(3,1) * t991 + rSges(3,2) * t931) * t858 + (rSges(3,1) * t931 + rSges(3,2) * t991) * t849) * m(3)) * t1015 + (t881 + 0.2e1 * t889 + t1013) * t702 + t729 * t696 - t753 * t928;
t924 = t766 * t954;
t697 = (-t877 * t987 + (pkin(1) * t727 + (0.2e1 * t712 * t949 - t924) * pkin(2)) * pkin(1)) * t727 * t903 + (-pkin(2) * t987 + (-t940 * t949 - t924) * pkin(1)) * t721 * t919 + ((pkin(1) * t940 * t954 + pkin(2) * t766) * t861 + pkin(1) * t948) * t836 * t766 * t957;
t703 = (-t724 * t993 + (-t861 * t1018 - t981) * pkin(2)) * t957;
t882 = t894 + t1012;
t888 = ((-t867 + (-rSges(3,1) * t861 + rSges(3,2) * t852) * m(3)) * t862 + t791 * t853) * pkin(1);
t730 = t888 + t882;
t754 = -t852 * Icges(3,5) - Icges(3,6) * t861 + (rSges(3,1) * t852 + rSges(3,2) * t861) * m(3) * (pkin(1) * t853 + rSges(3,3));
t885 = (-t795 * t852 - t796 * t861) * t763;
t891 = -t791 * t862 - t853 * t867;
t927 = t852 * t981;
t930 = t940 * t948;
t943 = t885 - (t852 * t861 * t1011 + (-0.4e1 * t834 + 0.2e1) * t797) * t978 + (t891 * t984 + ((-rSges(3,1) * t990 + rSges(3,2) * t930) * t861 + (rSges(3,1) * t930 + rSges(3,2) * t990) * t852) * m(3)) * t1015 + (t881 + 0.2e1 * t888 + t1012) * t703 + t730 * t697 - t754 * t927;
t770 = -t795 * t855 + t796 * t846;
t838 = -qJ(3,3) + qJ(2,3);
t799 = sin(t838);
t808 = cos(t838);
t913 = t825 * (t884 * t695 + t728 * t701 - t770 * t929 + t887 - (-0.2e1 * t971 - t974) * t980 + (((-t808 / 0.2e1 + t1004) * rSges(3,2) + (t799 / 0.2e1 + t1007) * rSges(3,1)) * m(3) - t893) * t998);
t771 = -t795 * t858 + t796 * t849;
t840 = qJ(2,2) - qJ(3,2);
t802 = sin(t840);
t810 = cos(t840);
t912 = t826 * (t883 * t696 + t729 * t702 - t771 * t928 + t886 - (-0.2e1 * t970 - t973) * t979 + (((-t810 / 0.2e1 + t1003) * rSges(3,2) + (t802 / 0.2e1 + t1006) * rSges(3,1)) * m(3) - t892) * t997);
t772 = -t795 * t861 + t796 * t852;
t842 = -qJ(3,1) + qJ(2,1);
t805 = sin(t842);
t812 = cos(t842);
t911 = t827 * (t882 * t697 + t730 * t703 - t772 * t927 + t885 - (-0.2e1 * t969 - t972) * t978 + (((-t812 / 0.2e1 + t1002) * rSges(3,2) + (t805 / 0.2e1 + t1005) * rSges(3,1)) * m(3) - t891) * t996);
t790 = -t939 * m(3) - Icges(3,3);
t910 = (t770 * t695 + t752 * t701 - t790 * t929 + (t974 / 0.2e1 + t971) * t1016 + ((t808 / 0.2e1 + t1004) * rSges(3,2) + (-t799 / 0.2e1 + t1007) * rSges(3,1)) * m(3) * t998) * t829;
t909 = (t771 * t696 + t753 * t702 - t790 * t928 + (t973 / 0.2e1 + t970) * t1017 + ((t810 / 0.2e1 + t1003) * rSges(3,2) + (-t802 / 0.2e1 + t1006) * rSges(3,1)) * m(3) * t997) * t832;
t908 = (t772 * t697 + t754 * t703 - t790 * t927 + (t972 / 0.2e1 + t969) * t1018 + ((t812 / 0.2e1 + t1002) * rSges(3,2) + (-t805 / 0.2e1 + t1005) * rSges(3,1)) * m(3) * t996) * t835;
t900 = t945 * t829 * t825;
t899 = t944 * t832 * t826;
t898 = t943 * t835 * t827;
t897 = t830 * t913;
t896 = t833 * t912;
t895 = t836 * t911;
t1 = [(t756 * t900 + t758 * t899 + t760 * t898) * t880 + (t815 * t908 + t814 * t909 + t813 * t910 + (t740 * t897 + t741 * t896 + t742 * t895) * t880) * t878; (t755 * t900 + t757 * t899 + t759 * t898) * t880 + (t818 * t908 + t817 * t909 + t816 * t910 + (t737 * t897 + t738 * t896 + t739 * t895) * t880) * t878; (-t943 * t966 - t944 * t967 - t945 * t968 + (t829 * t767 * t913 + t832 * t768 * t912 + t835 * t769 * t911) * t878) * t880;];
taucX  = t1;
