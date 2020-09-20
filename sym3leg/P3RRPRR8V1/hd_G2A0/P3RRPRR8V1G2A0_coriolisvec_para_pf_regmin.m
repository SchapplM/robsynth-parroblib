% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x13]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:34
% EndTime: 2020-08-06 19:59:38
% DurationCPUTime: 4.82s
% Computational Cost: add. (19006->329), mult. (42918->666), div. (2031->24), fcn. (28644->35), ass. (0->316)
t851 = pkin(4) + qJ(3,3);
t834 = 0.1e1 / t851;
t887 = (t851 ^ 2);
t835 = 0.1e1 / t887;
t836 = t834 * t835;
t852 = pkin(4) + qJ(3,2);
t837 = 0.1e1 / t852;
t889 = (t852 ^ 2);
t838 = 0.1e1 / t889;
t839 = t837 * t838;
t853 = pkin(4) + qJ(3,1);
t840 = 0.1e1 / t853;
t891 = (t853 ^ 2);
t841 = 0.1e1 / t891;
t842 = t840 * t841;
t1077 = 2 * pkin(1);
t850 = cos(pkin(5));
t1068 = pkin(2) * t850;
t807 = pkin(1) + t1068;
t871 = xDP(2);
t1013 = t871 * t807;
t864 = cos(qJ(1,3));
t1019 = t851 * t864;
t863 = cos(qJ(2,3));
t1026 = t807 * t863;
t1067 = pkin(2) * t871;
t849 = sin(pkin(5));
t803 = t849 * t1067;
t872 = xDP(1);
t1066 = pkin(2) * t872;
t804 = t849 * t1066;
t858 = sin(qJ(1,3));
t808 = t858 * t851;
t854 = legFrame(3,2);
t820 = sin(t854);
t823 = cos(t854);
t857 = sin(qJ(2,3));
t870 = xDP(3);
t937 = t807 * t872;
t1069 = pkin(2) * t849;
t1000 = t870 * t1069;
t953 = t864 * t1000;
t958 = t858 * t804;
t959 = t858 * t803;
t739 = ((t858 * t937 + t803) * t863 - t872 * t1019) * t823 + ((-t858 * t1013 + t804) * t863 + t871 * t1019) * t820 + t870 * (t864 * t1026 + t808) + ((-t958 + t1013) * t823 + (t959 + t937) * t820 - t953) * t857;
t1086 = 0.2e1 * t739;
t866 = cos(qJ(1,2));
t1018 = t852 * t866;
t865 = cos(qJ(2,2));
t1025 = t807 * t865;
t860 = sin(qJ(1,2));
t809 = t860 * t852;
t855 = legFrame(2,2);
t821 = sin(t855);
t824 = cos(t855);
t859 = sin(qJ(2,2));
t952 = t866 * t1000;
t956 = t860 * t804;
t957 = t860 * t803;
t740 = ((t860 * t937 + t803) * t865 - t872 * t1018) * t824 + ((-t860 * t1013 + t804) * t865 + t871 * t1018) * t821 + t870 * (t866 * t1025 + t809) + ((-t956 + t1013) * t824 + (t957 + t937) * t821 - t952) * t859;
t1085 = 0.2e1 * t740;
t868 = cos(qJ(1,1));
t1017 = t853 * t868;
t867 = cos(qJ(2,1));
t1024 = t807 * t867;
t862 = sin(qJ(1,1));
t810 = t862 * t853;
t856 = legFrame(1,2);
t822 = sin(t856);
t825 = cos(t856);
t861 = sin(qJ(2,1));
t951 = t868 * t1000;
t954 = t862 * t804;
t955 = t862 * t803;
t741 = ((t862 * t937 + t803) * t867 - t872 * t1017) * t825 + ((-t862 * t1013 + t804) * t867 + t871 * t1017) * t822 + t870 * (t868 * t1024 + t810) + ((-t954 + t1013) * t825 + (t955 + t937) * t822 - t951) * t861;
t1084 = 0.2e1 * t741;
t1005 = pkin(1) * t872 + t850 * t1066;
t1006 = pkin(1) * t871 + t850 * t1067;
t1023 = t807 * t870;
t745 = (t864 * t1023 + t803 * t823 + t804 * t820 + (-t1013 * t820 + t823 * t937) * t858) * t863 + t857 * ((-t958 + t1006) * t823 + (t959 + t1005) * t820 - t953);
t781 = t820 * t872 + t823 * t871;
t1083 = t745 * t781;
t746 = (t866 * t1023 + t803 * t824 + t804 * t821 + (-t1013 * t821 + t824 * t937) * t860) * t865 + t859 * ((-t956 + t1006) * t824 + (t957 + t1005) * t821 - t952);
t782 = t821 * t872 + t824 * t871;
t1082 = t746 * t782;
t747 = (t868 * t1023 + t803 * t825 + t804 * t822 + (-t1013 * t822 + t825 * t937) * t862) * t867 + t861 * ((-t954 + t1006) * t825 + (t955 + t1005) * t822 - t951);
t783 = t822 * t872 + t825 * t871;
t1081 = t747 * t783;
t826 = t863 * pkin(1);
t831 = qJ(2,3) + pkin(5);
t1080 = t826 + pkin(2) * cos(t831);
t827 = t865 * pkin(1);
t832 = qJ(2,2) + pkin(5);
t1079 = t827 + pkin(2) * cos(t832);
t828 = t867 * pkin(1);
t833 = qJ(2,1) + pkin(5);
t1078 = t828 + pkin(2) * cos(t833);
t784 = 0.1e1 / t1080;
t788 = 0.1e1 / t1079;
t792 = 0.1e1 / t1078;
t878 = t1080 ^ 2;
t785 = 0.1e1 / t878;
t881 = t1079 ^ 2;
t789 = 0.1e1 / t881;
t884 = t1078 ^ 2;
t793 = 0.1e1 / t884;
t877 = pkin(1) ^ 2;
t1076 = -2 * t877;
t846 = t863 ^ 2;
t817 = 0.2e1 * t846 - 0.1e1;
t847 = t865 ^ 2;
t818 = 0.2e1 * t847 - 0.1e1;
t848 = t867 ^ 2;
t819 = 0.2e1 * t848 - 0.1e1;
t778 = t781 ^ 2;
t1075 = pkin(1) * t778;
t779 = t782 ^ 2;
t1074 = pkin(1) * t779;
t780 = t783 ^ 2;
t1073 = pkin(1) * t780;
t1022 = t849 * t857;
t775 = -pkin(2) * t1022 + t1026;
t766 = 0.1e1 / t775;
t1053 = t745 * t766;
t1059 = t739 * t836;
t876 = pkin(2) ^ 2;
t1004 = -t876 - t877;
t802 = t1068 * t1077 - t1004;
t987 = t863 * t1053;
t988 = t857 * t1053;
t989 = t835 * t1053;
t718 = -t1053 * t1059 + (-(t988 * t1069 - t807 * t987 + t739) * t989 + t778 * t785 * t802 / (t826 + (t850 * t863 - t1022) * pkin(2))) * t834;
t1065 = t718 * qJ(3,3);
t1064 = t718 * t834;
t1021 = t849 * t859;
t776 = -pkin(2) * t1021 + t1025;
t768 = 0.1e1 / t776;
t1052 = t746 * t768;
t1058 = t740 * t839;
t984 = t865 * t1052;
t985 = t859 * t1052;
t986 = t838 * t1052;
t719 = -t1052 * t1058 + (-(t985 * t1069 - t807 * t984 + t740) * t986 + t779 * t789 * t802 / (t827 + (t850 * t865 - t1021) * pkin(2))) * t837;
t1063 = t719 * qJ(3,2);
t1062 = t719 * t837;
t1020 = t849 * t861;
t777 = -pkin(2) * t1020 + t1024;
t770 = 0.1e1 / t777;
t1051 = t747 * t770;
t1057 = t741 * t842;
t981 = t867 * t1051;
t982 = t861 * t1051;
t983 = t841 * t1051;
t720 = -t1051 * t1057 + (-(t982 * t1069 - t807 * t981 + t741) * t983 + t780 * t793 * t802 / (t828 + (t850 * t867 - t1020) * pkin(2))) * t840;
t1061 = t720 * qJ(3,1);
t1060 = t720 * t840;
t767 = 0.1e1 / t775 ^ 2;
t1056 = t745 ^ 2 * t767;
t769 = 0.1e1 / t776 ^ 2;
t1055 = t746 ^ 2 * t769;
t771 = 0.1e1 / t777 ^ 2;
t1054 = t747 ^ 2 * t771;
t1050 = t778 * t834;
t1049 = t779 * t837;
t1048 = t780 * t840;
t1047 = t781 * t784;
t1046 = t782 * t788;
t1045 = t783 * t792;
t1044 = t784 * t820;
t1043 = t784 * t823;
t1042 = t784 * t835;
t1041 = t785 * t863;
t796 = pkin(1) * t857 + pkin(2) * sin(t831);
t1040 = t784 * t785 * t796;
t1039 = t788 * t821;
t1038 = t788 * t824;
t1037 = t788 * t838;
t1036 = t789 * t865;
t797 = pkin(1) * t859 + pkin(2) * sin(t832);
t1035 = t788 * t789 * t797;
t1034 = t792 * t822;
t1033 = t792 * t825;
t1032 = t792 * t841;
t1031 = t793 * t867;
t798 = pkin(1) * t861 + pkin(2) * sin(t833);
t1030 = t792 * t793 * t798;
t1029 = t807 * t858;
t1028 = t807 * t860;
t1027 = t807 * t862;
t1016 = t857 * t863;
t1015 = t859 * t865;
t1014 = t861 * t867;
t730 = -t1041 * t1075 + t989 * t1086;
t873 = 0.2e1 * qJ(2,3);
t999 = pkin(2) * t1077;
t895 = -(-t802 * t1047 + t796 * t1053) * t834 * t1047 + (-t796 * t851 * t1047 * t835 - t1080 * t1059 - ((-t876 * cos(0.2e1 * t831) - cos(t873) * t877 - (2 * t887) + (-cos(t873 + pkin(5)) - t850) * t999 + t1004) * t1053 + t1080 * t1086) * t836 / 0.2e1) * t1053;
t950 = t834 * t988;
t924 = t950 * t1047;
t977 = t857 * t1040;
t947 = t778 * t977;
t1012 = (qJ(3,3) ^ 2 + t846 * t877) * t718 + (-qJ(3,3) * t947 - t863 * t895) * pkin(1) + t863 * t924 * t1076 + t730 * qJ(3,3);
t731 = -t1036 * t1074 + t986 * t1085;
t874 = 0.2e1 * qJ(2,2);
t894 = -(-t802 * t1046 + t797 * t1052) * t837 * t1046 + (-t797 * t852 * t1046 * t838 - t1079 * t1058 - ((-t876 * cos(0.2e1 * t832) - cos(t874) * t877 - (2 * t889) + (-cos(pkin(5) + t874) - t850) * t999 + t1004) * t1052 + t1079 * t1085) * t839 / 0.2e1) * t1052;
t949 = t837 * t985;
t923 = t949 * t1046;
t976 = t859 * t1035;
t946 = t779 * t976;
t1011 = (qJ(3,2) ^ 2 + t847 * t877) * t719 + (-qJ(3,2) * t946 - t865 * t894) * pkin(1) + t865 * t923 * t1076 + t731 * qJ(3,2);
t732 = -t1031 * t1073 + t983 * t1084;
t875 = 0.2e1 * qJ(2,1);
t893 = -(-t802 * t1045 + t798 * t1051) * t840 * t1045 + (-t798 * t853 * t1045 * t841 - t1078 * t1057 - ((-t876 * cos(0.2e1 * t833) - cos(t875) * t877 - (2 * t891) + (-cos(pkin(5) + t875) - t850) * t999 + t1004) * t1051 + t1078 * t1084) * t842 / 0.2e1) * t1051;
t948 = t840 * t982;
t922 = t948 * t1045;
t975 = t861 * t1030;
t945 = t780 * t975;
t1010 = (qJ(3,1) ^ 2 + t848 * t877) * t720 + (-qJ(3,1) * t945 - t867 * t893) * pkin(1) + t867 * t922 * t1076 + t732 * qJ(3,1);
t995 = t835 * t1056;
t1009 = -qJ(3,3) * t995 + t924 * t1077 - t718 * t826 + t895;
t993 = t838 * t1055;
t1008 = -qJ(3,2) * t993 + t923 * t1077 - t719 * t827 + t894;
t991 = t841 * t1054;
t1007 = -qJ(3,1) * t991 + t922 * t1077 - t720 * t828 + t893;
t1003 = t858 * t1069;
t1002 = t860 * t1069;
t1001 = t862 * t1069;
t998 = t766 * t1064;
t997 = t768 * t1062;
t996 = t770 * t1060;
t994 = t836 * t1056;
t992 = t839 * t1055;
t990 = t842 * t1054;
t980 = t778 / t878 ^ 2 * t796;
t979 = t779 / t881 ^ 2 * t797;
t978 = t780 / t884 ^ 2 * t798;
t974 = t718 * t1044;
t973 = t719 * t1039;
t972 = t720 * t1034;
t971 = t718 * t1043;
t970 = t719 * t1038;
t969 = t720 * t1033;
t968 = t857 ^ 2 * t1064;
t967 = t859 ^ 2 * t1062;
t966 = t861 ^ 2 * t1060;
t965 = t1012 * t766;
t964 = t1011 * t768;
t963 = t1010 * t770;
t962 = (-pkin(1) * t947 + 0.2e1 * t1065 + t730) * t834;
t961 = (-pkin(1) * t946 + 0.2e1 * t1063 + t731) * t837;
t960 = (-pkin(1) * t945 + 0.2e1 * t1061 + t732) * t840;
t944 = t1016 * t1042;
t943 = t1015 * t1037;
t942 = t1014 * t1032;
t941 = t1016 * t1064;
t940 = t1015 * t1062;
t939 = t1014 * t1060;
t936 = t766 * t962;
t935 = t768 * t961;
t934 = t770 * t960;
t930 = t817 * t1042 * t1083;
t929 = t818 * t1037 * t1082;
t928 = t819 * t1032 * t1081;
t927 = t944 * t1056;
t926 = t943 * t1055;
t925 = t942 * t1054;
t921 = 0.2e1 * t944 * t1083;
t920 = 0.2e1 * t943 * t1082;
t919 = 0.2e1 * t942 * t1081;
t918 = t784 * (t1040 * t1075 - t857 * t1065 + (pkin(1) * t987 - 0.2e1 * t739) * t834 * t950);
t917 = t788 * (t1035 * t1074 - t859 * t1063 + (pkin(1) * t984 - 0.2e1 * t740) * t837 * t949);
t916 = t792 * (t1030 * t1073 - t861 * t1061 + (pkin(1) * t981 - 0.2e1 * t741) * t840 * t948);
t915 = (t863 * t1040 - t785 * t857) * t1050;
t914 = (t977 + t1041) * t1050;
t913 = (t865 * t1035 - t789 * t859) * t1049;
t912 = (t976 + t1036) * t1049;
t911 = (t867 * t1030 - t793 * t861) * t1048;
t910 = (t975 + t1031) * t1048;
t909 = t766 * t915;
t908 = t766 * t914;
t907 = t768 * t913;
t906 = t768 * t912;
t905 = t770 * t911;
t904 = t770 * t910;
t903 = t766 * t968 + t767 * t921;
t902 = t768 * t967 + t769 * t920;
t901 = t770 * t966 + t771 * t919;
t900 = t857 * t974 + t859 * t973 + t861 * t972;
t899 = t857 * t971 + t859 * t970 + t861 * t969;
t898 = 0.2e1 * t766 * t941 + 0.2e1 * t767 * t930;
t897 = 0.2e1 * t768 * t940 + 0.2e1 * t769 * t929;
t896 = 0.2e1 * t770 * t939 + 0.2e1 * t771 * t928;
t774 = t867 * t1069 + t807 * t861;
t773 = t865 * t1069 + t807 * t859;
t772 = t863 * t1069 + t807 * t857;
t765 = t1078 * t868 + t810;
t764 = t1079 * t866 + t809;
t763 = t1080 * t864 + t808;
t762 = t777 * t862 - t1017;
t761 = t776 * t860 - t1018;
t760 = t775 * t858 - t1019;
t759 = (t825 * t1027 + t822 * t1069) * t867 + (-t825 * t1001 + t807 * t822) * t861;
t758 = (t824 * t1028 + t821 * t1069) * t865 + (-t824 * t1002 + t807 * t821) * t859;
t757 = (t823 * t1029 + t820 * t1069) * t863 + (-t823 * t1003 + t807 * t820) * t857;
t756 = (-t822 * t1027 + t825 * t1069) * t867 + t861 * (t822 * t1001 + t807 * t825);
t755 = (-t821 * t1028 + t824 * t1069) * t865 + t859 * (t821 * t1002 + t807 * t824);
t754 = (-t820 * t1029 + t823 * t1069) * t863 + t857 * (t820 * t1003 + t807 * t823);
t753 = -t762 * t822 + t774 * t825;
t752 = t762 * t825 + t774 * t822;
t751 = -t761 * t821 + t773 * t824;
t750 = t761 * t824 + t773 * t821;
t749 = -t760 * t820 + t772 * t823;
t748 = t760 * t823 + t772 * t820;
t735 = t819 * t991;
t734 = t818 * t993;
t733 = t817 * t995;
t1 = [t757 * t998 + t758 * t997 + t759 * t996, 0, 0, t903 * t757 + t902 * t758 + t901 * t759 - t820 * t927 - t821 * t926 - t822 * t925, -t735 * t1034 - t734 * t1039 - t733 * t1044 + t757 * t898 + t758 * t897 + t759 * t896, t757 * t908 + t758 * t906 + t759 * t904 + t900, t757 * t909 + t758 * t907 + t759 * t905 + t863 * t974 + t865 * t973 + t867 * t972, t820 * t980 + t821 * t979 + t822 * t978, 0, 0, -pkin(1) * t900 - t748 * t994 - t750 * t992 - t752 * t990 + t757 * t936 + t758 * t935 + t759 * t934, (t1007 * t752 + t759 * t963) * t840 + (t1008 * t750 + t758 * t964) * t837 + (t1009 * t748 + t757 * t965) * t834 + (t820 * t918 + t821 * t917 + t822 * t916) * pkin(1), 0; t754 * t998 + t755 * t997 + t756 * t996, 0, 0, t903 * t754 + t902 * t755 + t901 * t756 - t823 * t927 - t824 * t926 - t825 * t925, -t735 * t1033 - t734 * t1038 - t733 * t1043 + t754 * t898 + t755 * t897 + t756 * t896, t754 * t908 + t755 * t906 + t756 * t904 + t899, t754 * t909 + t755 * t907 + t756 * t905 + t863 * t971 + t865 * t970 + t867 * t969, t823 * t980 + t824 * t979 + t825 * t978, 0, 0, -pkin(1) * t899 - t749 * t994 - t751 * t992 - t753 * t990 + t754 * t936 + t755 * t935 + t756 * t934, (t1007 * t753 + t756 * t963) * t840 + (t1008 * t751 + t755 * t964) * t837 + (t1009 * t749 + t754 * t965) * t834 + (t823 * t918 + t824 * t917 + t825 * t916) * pkin(1), 0; t868 * t1060 + t866 * t1062 + t864 * t1064, 0, 0, (t770 * t919 + t966) * t868 + (t768 * t920 + t967) * t866 + (t766 * t921 + t968) * t864, 0.2e1 * (t770 * t928 + t939) * t868 + 0.2e1 * (t768 * t929 + t940) * t866 + 0.2e1 * (t766 * t930 + t941) * t864, t864 * t914 + t866 * t912 + t868 * t910, t864 * t915 + t866 * t913 + t868 * t911, 0, 0, 0, -t763 * t994 - t764 * t992 - t765 * t990 + t864 * t962 + t866 * t961 + t868 * t960, (t1007 * t765 + t1010 * t868) * t840 + (t1008 * t764 + t1011 * t866) * t837 + (t1009 * t763 + t1012 * t864) * t834, 0;];
tau_reg  = t1;
