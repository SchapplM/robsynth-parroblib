% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR8V2G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:02:31
% EndTime: 2020-08-06 21:02:38
% DurationCPUTime: 6.47s
% Computational Cost: add. (27795->456), mult. (34215->719), div. (5334->15), fcn. (25467->62), ass. (0->340)
t891 = rSges(2,3) + pkin(5);
t864 = qJ(2,1) + pkin(7);
t832 = cos(t864);
t1052 = pkin(3) * t832;
t889 = cos(qJ(2,1));
t846 = t889 * pkin(2);
t810 = t846 + pkin(1);
t1047 = pkin(5) + qJ(3,1);
t853 = -pkin(6) - t1047;
t884 = sin(qJ(1,1));
t890 = cos(qJ(1,1));
t756 = t810 * t884 + t853 * t890;
t759 = t810 * t890 - t853 * t884;
t878 = legFrame(1,3);
t836 = sin(t878);
t839 = cos(t878);
t765 = t836 * t890 + t839 * t884;
t734 = t1052 * t765 + t756 * t839 + t759 * t836;
t897 = xDP(2);
t1036 = t734 * t897;
t762 = -t836 * t884 + t839 * t890;
t731 = t1052 * t762 - t756 * t836 + t759 * t839;
t898 = xDP(1);
t1039 = t731 * t898;
t842 = 0.1e1 / t853;
t784 = t846 + t1052;
t775 = 0.1e1 / t784;
t896 = xDP(3);
t1027 = t775 * t896;
t1095 = 0.2e1 * pkin(1);
t883 = sin(qJ(2,1));
t1056 = pkin(2) * t883;
t821 = sin(t864);
t781 = pkin(3) * t821 + t1056;
t806 = 0.2e1 * t864;
t797 = sin(t806);
t912 = pkin(3) ^ 2;
t1002 = t912 * t797;
t906 = 0.2e1 * qJ(2,1);
t863 = t906 + pkin(7);
t820 = sin(t863);
t998 = 0.2e1 * pkin(2) * pkin(3);
t868 = sin(t906);
t914 = pkin(2) ^ 2;
t999 = t914 * t868;
t924 = t820 * t998 + t1002 + t999;
t743 = t1095 * t781 + t924;
t964 = t743 * t1027;
t713 = (t964 / 0.2e1 + t1036 + t1039) * t842;
t862 = qJ(2,2) + pkin(7);
t828 = cos(t862);
t1053 = pkin(3) * t828;
t887 = cos(qJ(2,2));
t845 = t887 * pkin(2);
t809 = t845 + pkin(1);
t1046 = pkin(5) + qJ(3,2);
t852 = -pkin(6) - t1046;
t882 = sin(qJ(1,2));
t888 = cos(qJ(1,2));
t755 = t809 * t882 + t852 * t888;
t758 = t809 * t888 - t852 * t882;
t877 = legFrame(2,3);
t835 = sin(t877);
t838 = cos(t877);
t764 = t835 * t888 + t838 * t882;
t733 = t1053 * t764 + t755 * t838 + t758 * t835;
t1037 = t733 * t897;
t761 = -t835 * t882 + t838 * t888;
t730 = t1053 * t761 - t755 * t835 + t758 * t838;
t1040 = t730 * t898;
t841 = 0.1e1 / t852;
t783 = t845 + t1053;
t772 = 0.1e1 / t783;
t1029 = t772 * t896;
t881 = sin(qJ(2,2));
t1057 = pkin(2) * t881;
t819 = sin(t862);
t780 = pkin(3) * t819 + t1057;
t904 = 0.2e1 * qJ(2,2);
t867 = sin(t904);
t1000 = t914 * t867;
t805 = 0.2e1 * t862;
t796 = sin(t805);
t1003 = t912 * t796;
t861 = t904 + pkin(7);
t818 = sin(t861);
t925 = t818 * t998 + t1000 + t1003;
t742 = t1095 * t780 + t925;
t965 = t742 * t1029;
t712 = (t965 / 0.2e1 + t1037 + t1040) * t841;
t860 = qJ(2,3) + pkin(7);
t825 = cos(t860);
t1054 = pkin(3) * t825;
t885 = cos(qJ(2,3));
t844 = t885 * pkin(2);
t808 = t844 + pkin(1);
t1045 = pkin(5) + qJ(3,3);
t851 = -pkin(6) - t1045;
t880 = sin(qJ(1,3));
t886 = cos(qJ(1,3));
t754 = t808 * t880 + t851 * t886;
t757 = t808 * t886 - t851 * t880;
t876 = legFrame(3,3);
t834 = sin(t876);
t837 = cos(t876);
t763 = t834 * t886 + t837 * t880;
t732 = t1054 * t763 + t754 * t837 + t757 * t834;
t1038 = t732 * t897;
t760 = -t834 * t880 + t837 * t886;
t729 = t1054 * t760 - t754 * t834 + t757 * t837;
t1041 = t729 * t898;
t840 = 0.1e1 / t851;
t782 = t844 + t1054;
t769 = 0.1e1 / t782;
t1031 = t769 * t896;
t879 = sin(qJ(2,3));
t1058 = pkin(2) * t879;
t817 = sin(t860);
t779 = pkin(3) * t817 + t1058;
t902 = 0.2e1 * qJ(2,3);
t866 = sin(t902);
t1001 = t914 * t866;
t804 = 0.2e1 * t860;
t795 = sin(t804);
t1004 = t912 * t795;
t859 = t902 + pkin(7);
t816 = sin(t859);
t926 = t816 * t998 + t1001 + t1004;
t741 = t1095 * t779 + t926;
t966 = t741 * t1031;
t711 = (t966 / 0.2e1 + t1038 + t1041) * t840;
t1094 = -t896 / 0.2e1;
t874 = cos(pkin(7));
t803 = pkin(3) * t874 + pkin(2);
t770 = 0.1e1 / t782 ^ 2;
t872 = t896 ^ 2;
t1093 = t770 * t872;
t773 = 0.1e1 / t783 ^ 2;
t1092 = t773 * t872;
t776 = 0.1e1 / t784 ^ 2;
t1091 = t776 * t872;
t1078 = m(2) * rSges(2,2);
t812 = rSges(2,1) * t1078 - Icges(2,4);
t869 = cos(t902);
t1013 = t812 * t869;
t798 = cos(t804);
t1077 = m(3) * rSges(3,2);
t811 = -rSges(3,1) * t1077 + Icges(3,4);
t1016 = t811 * t798;
t907 = rSges(3,2) ^ 2;
t909 = rSges(3,1) ^ 2;
t787 = m(3) * (-t907 + t909) - Icges(3,1) + Icges(3,2);
t1022 = t787 * t795;
t908 = rSges(2,2) ^ 2;
t910 = rSges(2,1) ^ 2;
t778 = m(3) * t914 + (-t908 + t910) * m(2) + Icges(2,2) - Icges(2,1);
t1026 = t778 * t866;
t873 = sin(pkin(7));
t1010 = t873 * t879;
t973 = pkin(3) * t1010;
t1088 = t803 * t885 - t973;
t748 = 0.1e1 / t1088;
t1051 = pkin(3) * t873;
t751 = t1051 * t885 + t803 * t879;
t723 = (t748 * t751 * t896 + t760 * t898 + t763 * t897) * t840;
t848 = rSges(3,3) + t1045;
t1044 = t723 * t848;
t1061 = pkin(2) * t816;
t1064 = pkin(1) * t825;
t1070 = t812 / 0.2e1;
t1071 = -t811 / 0.2e1;
t1072 = -t787 / 0.2e1;
t1073 = -t778 / 0.2e1;
t1086 = -0.2e1 * t711;
t1069 = -0.3e1 / 0.4e1 * t914;
t1076 = -t723 / 0.4e1;
t1055 = pkin(2) * t912;
t913 = pkin(2) * t914;
t1081 = -0.2e1 * t913 - 0.4e1 * t1055;
t1048 = t874 * pkin(2);
t814 = pkin(3) * t1048;
t993 = t814 + t914 / 0.2e1;
t1082 = -0.4e1 * pkin(1) * (t912 / 0.2e1 + t993);
t720 = pkin(1) * t723;
t702 = -t720 - (-t1041 / 0.2e1 - t1038 / 0.2e1 - t966 / 0.4e1) * t840;
t915 = pkin(1) ^ 2;
t717 = (t851 ^ 2 + t915) * t723;
t738 = t851 * t1031;
t900 = 0.2e1 * pkin(7);
t822 = cos(t900 + qJ(2,3));
t826 = cos(qJ(2,3) - pkin(7));
t847 = t912 + t914;
t901 = 0.3e1 * qJ(2,3);
t911 = pkin(3) * t912;
t951 = 0.3e1 / 0.4e1 * t914;
t952 = 0.3e1 / 0.4e1 * t912;
t1049 = pkin(3) * t914;
t953 = -0.2e1 * t911 - 0.4e1 * t1049;
t1032 = t769 * t840;
t956 = -t1032 / 0.2e1;
t982 = -0.2e1 * t1049;
t983 = -0.2e1 * t1055;
t988 = -0.6e1 * t914 - 0.3e1 * t912;
t824 = cos(t859);
t992 = t824 + t874;
t693 = (t1081 * t885 + t822 * t983 + t825 * t953 + t826 * t982 + t1082) * t956 * t1093 + (t926 * t770 * t851 * t1094 + ((-t911 * cos(0.3e1 * t860) - t913 * cos(t901)) * t1076 - (t1004 / 0.2e1 + t1001 / 0.2e1) * t738 + (-(-cos(t901 + pkin(7)) - t826) * t723 * t951 + (pkin(1) * t1086 + t1076 * t988 + t717) * t825) * pkin(3) + ((-t1069 * t723 + t717) * t885 - (-cos(t900 + t901) - t822 - 0.2e1 * t885) * t723 * t952 + (-0.2e1 * t702 * t992 - t738 * t816) * pkin(3) - (pkin(3) * t992 + t1095 * t885) * t711) * pkin(2) + (-t711 / 0.2e1 - t702) * (t798 * t912 + t869 * t914 + t847)) * t769) * t840 * t723;
t1007 = t912 * t723;
t1023 = (0.2e1 * t814 + t847) * t872;
t1035 = t748 * t840;
t858 = t874 ^ 2;
t1083 = -0.2e1 * (t858 - 0.1e1 / 0.2e1) * t912 - 0.2e1 * t993;
t705 = -t720 + t711;
t696 = t770 * t840 / (t844 + (t874 * t885 - t1010) * pkin(3)) * t1023 - (t858 * t1007 - t705 * t973 - t1007 + t1088 * t711 + (t723 * t1083 * t885 - t803 * (-0.2e1 * t723 * t973 - t705)) * t885) * t723 * t1035;
t1067 = m(3) * t848;
t789 = rSges(3,2) * t1067 - Icges(3,6);
t792 = rSges(3,1) * t1067 - Icges(3,5);
t1079 = m(2) * rSges(2,1);
t938 = -t1079 * t891 + Icges(2,5);
t957 = t891 * t1078 - Icges(2,6);
t979 = pkin(2) * t1067;
t726 = -(t789 * t873 - t792 * t874 + t938 - t979) * t879 + (t789 * t874 + t792 * t873 + t957) * t885;
t745 = -rSges(3,1) * t825 + rSges(3,2) * t817 - t808;
t1068 = m(2) * t891;
t801 = rSges(2,2) * t1068 - Icges(2,6);
t802 = t873 * pkin(2) * t1077;
t899 = 0.2e1 * t915;
t989 = t908 + t910;
t927 = -((2 * pkin(5) ^ 2) + t899 + ((4 * pkin(5) + 2 * rSges(2,3)) * rSges(2,3)) + t989) * m(2) / 0.2e1 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t802 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t928 = -t899 / 0.2e1 - t907 / 0.2e1 - t909 / 0.2e1 - t914 / 0.2e1;
t1080 = pkin(2) * m(3);
t807 = t1079 + t1080;
t934 = t1078 * t885 + t807 * t879;
t960 = t769 * t779 * t1093;
t970 = t1078 * t1095;
t980 = 0.2e1 * rSges(3,2);
t981 = 0.2e1 * rSges(3,1);
t984 = -0.2e1 * pkin(1) * t807;
t986 = rSges(3,2) * t1095;
t997 = rSges(2,1) * t1068 - Icges(2,5);
t996 = (t798 * t1072 + t869 * t1073 + t879 * t970 + t885 * t984 + t927) * t696 - t726 * t960 + 0.2e1 * (t866 * t1070 + t795 * t1071) * t696 + (-t745 * t693 + (rSges(3,2) * t1061 + t817 * t986 - t848 ^ 2 + t928 + (-pkin(2) * t992 - 0.2e1 * t1064) * rSges(3,1)) * t696) * m(3) - m(3) * t1044 * t1086 + ((-t792 * t825 + t789 * t817 - (t979 + t997) * t885 + t801 * t879) * t1031 - (-t1026 - t1022 + 0.2e1 * t1016 - 0.2e1 * t1013 - t934 * t1095 + ((-pkin(2) * t824 - t1064) * t980 + (-pkin(1) * t817 - t1061) * t981) * m(3)) * t723) * t1031;
t870 = cos(t904);
t1012 = t812 * t870;
t799 = cos(t805);
t1015 = t811 * t799;
t1021 = t787 * t796;
t1025 = t778 * t867;
t1009 = t873 * t881;
t972 = pkin(3) * t1009;
t1089 = t803 * t887 - t972;
t749 = 0.1e1 / t1089;
t752 = t1051 * t887 + t803 * t881;
t724 = (t749 * t752 * t896 + t761 * t898 + t764 * t897) * t841;
t849 = rSges(3,3) + t1046;
t1043 = t724 * t849;
t1060 = pkin(2) * t818;
t1063 = pkin(1) * t828;
t1085 = -0.2e1 * t712;
t1075 = -t724 / 0.4e1;
t721 = pkin(1) * t724;
t703 = -t721 - (-t1040 / 0.2e1 - t1037 / 0.2e1 - t965 / 0.4e1) * t841;
t718 = (t852 ^ 2 + t915) * t724;
t739 = t852 * t1029;
t823 = cos(t900 + qJ(2,2));
t829 = cos(qJ(2,2) - pkin(7));
t903 = 0.3e1 * qJ(2,2);
t1030 = t772 * t841;
t955 = -t1030 / 0.2e1;
t827 = cos(t861);
t991 = t827 + t874;
t694 = (t1081 * t887 + t823 * t983 + t828 * t953 + t829 * t982 + t1082) * t955 * t1092 + (t925 * t773 * t852 * t1094 + ((-t911 * cos(0.3e1 * t862) - t913 * cos(t903)) * t1075 - (t1003 / 0.2e1 + t1000 / 0.2e1) * t739 + (-(-cos(t903 + pkin(7)) - t829) * t724 * t951 + (pkin(1) * t1085 + t1075 * t988 + t718) * t828) * pkin(3) + ((-t1069 * t724 + t718) * t887 - (-cos(t900 + t903) - t823 - 0.2e1 * t887) * t724 * t952 + (-0.2e1 * t703 * t991 - t739 * t818) * pkin(3) - (pkin(3) * t991 + t1095 * t887) * t712) * pkin(2) + (-t712 / 0.2e1 - t703) * (t912 * t799 + t914 * t870 + t847)) * t772) * t841 * t724;
t1006 = t912 * t724;
t1034 = t749 * t841;
t706 = -t721 + t712;
t697 = t773 * t841 / (t845 + (t874 * t887 - t1009) * pkin(3)) * t1023 - (t858 * t1006 - t706 * t972 - t1006 + t1089 * t712 + (t724 * t1083 * t887 - t803 * (-0.2e1 * t724 * t972 - t706)) * t887) * t724 * t1034;
t1066 = m(3) * t849;
t790 = rSges(3,2) * t1066 - Icges(3,6);
t793 = rSges(3,1) * t1066 - Icges(3,5);
t978 = pkin(2) * t1066;
t727 = -(t790 * t873 - t793 * t874 + t938 - t978) * t881 + (t790 * t874 + t793 * t873 + t957) * t887;
t746 = -rSges(3,1) * t828 + rSges(3,2) * t819 - t809;
t933 = t1078 * t887 + t807 * t881;
t959 = t772 * t780 * t1092;
t995 = (t799 * t1072 + t870 * t1073 + t881 * t970 + t887 * t984 + t927) * t697 - t727 * t959 + 0.2e1 * (t867 * t1070 + t796 * t1071) * t697 + (-t746 * t694 + (rSges(3,2) * t1060 + t819 * t986 - t849 ^ 2 + t928 + (-pkin(2) * t991 - 0.2e1 * t1063) * rSges(3,1)) * t697) * m(3) - m(3) * t1043 * t1085 + ((-t793 * t828 + t790 * t819 - (t978 + t997) * t887 + t801 * t881) * t1029 - (-t1025 - t1021 + 0.2e1 * t1015 - 0.2e1 * t1012 - t933 * t1095 + ((-pkin(2) * t827 - t1063) * t980 + (-pkin(1) * t819 - t1060) * t981) * m(3)) * t724) * t1029;
t871 = cos(t906);
t1011 = t812 * t871;
t800 = cos(t806);
t1014 = t811 * t800;
t1020 = t787 * t797;
t1024 = t778 * t868;
t1008 = t873 * t883;
t971 = pkin(3) * t1008;
t1090 = t803 * t889 - t971;
t750 = 0.1e1 / t1090;
t753 = t1051 * t889 + t803 * t883;
t725 = (t750 * t753 * t896 + t762 * t898 + t765 * t897) * t842;
t850 = rSges(3,3) + t1047;
t1042 = t725 * t850;
t1059 = pkin(2) * t820;
t1062 = pkin(1) * t832;
t1084 = -0.2e1 * t713;
t1074 = -t725 / 0.4e1;
t722 = pkin(1) * t725;
t704 = -t722 - (-t1039 / 0.2e1 - t1036 / 0.2e1 - t964 / 0.4e1) * t842;
t719 = (t853 ^ 2 + t915) * t725;
t740 = t853 * t1027;
t831 = cos(qJ(2,1) + t900);
t833 = cos(qJ(2,1) - pkin(7));
t905 = 0.3e1 * qJ(2,1);
t1028 = t775 * t842;
t954 = -t1028 / 0.2e1;
t830 = cos(t863);
t990 = t830 + t874;
t695 = (t1081 * t889 + t831 * t983 + t832 * t953 + t833 * t982 + t1082) * t954 * t1091 + (t924 * t776 * t853 * t1094 + ((-t911 * cos(0.3e1 * t864) - t913 * cos(t905)) * t1074 - (t1002 / 0.2e1 + t999 / 0.2e1) * t740 + (-(-cos(t905 + pkin(7)) - t833) * t725 * t951 + (pkin(1) * t1084 + t1074 * t988 + t719) * t832) * pkin(3) + ((-t1069 * t725 + t719) * t889 - (-cos(t900 + t905) - t831 - 0.2e1 * t889) * t725 * t952 + (-0.2e1 * t704 * t990 - t740 * t820) * pkin(3) - (pkin(3) * t990 + t1095 * t889) * t713) * pkin(2) + (-t704 - t713 / 0.2e1) * (t912 * t800 + t914 * t871 + t847)) * t775) * t842 * t725;
t1005 = t912 * t725;
t1033 = t750 * t842;
t707 = -t722 + t713;
t698 = t776 * t842 / (t846 + (t874 * t889 - t1008) * pkin(3)) * t1023 - (t858 * t1005 - t707 * t971 - t1005 + t1090 * t713 + (t725 * t1083 * t889 - t803 * (-0.2e1 * t725 * t971 - t707)) * t889) * t725 * t1033;
t1065 = m(3) * t850;
t791 = rSges(3,2) * t1065 - Icges(3,6);
t794 = rSges(3,1) * t1065 - Icges(3,5);
t977 = pkin(2) * t1065;
t728 = -(t791 * t873 - t794 * t874 + t938 - t977) * t883 + (t791 * t874 + t794 * t873 + t957) * t889;
t747 = -rSges(3,1) * t832 + rSges(3,2) * t821 - t810;
t932 = t1078 * t889 + t807 * t883;
t958 = t775 * t781 * t1091;
t994 = (t800 * t1072 + t871 * t1073 + t883 * t970 + t889 * t984 + t927) * t698 - t728 * t958 + 0.2e1 * (t868 * t1070 + t797 * t1071) * t698 + (-t747 * t695 + (rSges(3,2) * t1059 + t821 * t986 - t850 ^ 2 + t928 + (-pkin(2) * t990 - 0.2e1 * t1062) * rSges(3,1)) * t698) * m(3) - m(3) * t1042 * t1084 + ((-t794 * t832 + t791 * t821 - (t977 + t997) * t889 + t801 * t883) * t1027 - (-t1024 - t1020 + 0.2e1 * t1014 - 0.2e1 * t1011 - t932 * t1095 + ((-pkin(2) * t830 - t1062) * t980 + (-pkin(1) * t821 - t1059) * t981) * m(3)) * t725) * t1027;
t944 = rSges(3,1) * t817 + rSges(3,2) * t825;
t976 = m(3) * (t1044 / 0.2e1 + (t944 + t1058) * t1031) * t723;
t943 = rSges(3,1) * t819 + rSges(3,2) * t828;
t975 = m(3) * (t1043 / 0.2e1 + (t943 + t1057) * t1029) * t724;
t942 = rSges(3,1) * t821 + rSges(3,2) * t832;
t974 = m(3) * (t1042 / 0.2e1 + (t942 + t1056) * t1027) * t725;
t690 = (-t696 * t745 - t693) * m(3);
t947 = t690 - 0.2e1 * t976;
t691 = (-t697 * t746 - t694) * m(3);
t946 = t691 - 0.2e1 * t975;
t692 = (-t698 * t747 - t695) * m(3);
t945 = t692 - 0.2e1 * t974;
t744 = 0.2e1 * t802 - t989 * m(2) - Icges(2,3) - Icges(3,3) + (-0.2e1 * rSges(3,1) * t1048 - t907 - t909 - t914) * m(3);
t1 = [-(t731 * t945 + t762 * t994) * t842 - (t730 * t946 + t761 * t995) * t841 - (t729 * t947 + t760 * t996) * t840; -(t734 * t945 + t765 * t994) * t842 - (t733 * t946 + t764 * t995) * t841 - (t732 * t947 + t763 * t996) * t840; -t996 * t751 * t1035 - t995 * t752 * t1034 - t994 * t753 * t1033 + (t698 * t728 - t744 * t958 - t725 * ((0.2e1 * t713 * t1056 + t942 * (-t722 - (-0.2e1 * t1036 - t964 - 0.2e1 * t1039) * t842)) * m(3) - (t1020 / 0.2e1 - t1014 + t1024 / 0.2e1 + t1011 + t932 * pkin(1) + (rSges(3,1) * t820 + rSges(3,2) * t830) * t1080) * t725)) * t775 + (t697 * t727 - t744 * t959 - t724 * ((0.2e1 * t712 * t1057 + t943 * (-t721 - (-0.2e1 * t1037 - t965 - 0.2e1 * t1040) * t841)) * m(3) - (t1021 / 0.2e1 - t1015 + t1025 / 0.2e1 + t1012 + t933 * pkin(1) + (rSges(3,1) * t818 + rSges(3,2) * t827) * t1080) * t724)) * t772 + (t696 * t726 - t744 * t960 - t723 * ((0.2e1 * t711 * t1058 + t944 * (-t720 - (-0.2e1 * t1038 - t966 - 0.2e1 * t1041) * t840)) * m(3) - (t1022 / 0.2e1 - t1016 + t1026 / 0.2e1 + t1013 + t934 * pkin(1) + (rSges(3,1) * t816 + rSges(3,2) * t824) * t1080) * t723)) * t769 + (t974 * t1028 + t692 * t954) * t743 + (t975 * t1030 + t691 * t955) * t742 + (t976 * t1032 + t690 * t956) * t741;];
taucX  = t1;
