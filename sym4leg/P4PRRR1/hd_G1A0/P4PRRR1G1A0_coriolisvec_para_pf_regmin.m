% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P4PRRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x11]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRR1G1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_regmin: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:14:44
% EndTime: 2020-03-02 20:14:55
% DurationCPUTime: 11.29s
% Computational Cost: add. (46867->291), mult. (40644->529), div. (5824->11), fcn. (35532->34), ass. (0->269)
t892 = pkin(7) + qJ(2,4);
t872 = qJ(3,4) + t892;
t858 = sin(t872);
t859 = cos(t872);
t870 = sin(t892);
t871 = cos(t892);
t806 = t858 * t871 - t859 * t870;
t1069 = 0.1e1 / t806;
t912 = xP(4);
t890 = sin(t912);
t891 = cos(t912);
t913 = koppelP(4,2);
t917 = koppelP(4,1);
t842 = t890 * t917 + t891 * t913;
t909 = xDP(4);
t911 = xDP(1);
t830 = t842 * t909 - t911;
t846 = -t890 * t913 + t891 * t917;
t910 = xDP(2);
t834 = t846 * t909 + t910;
t897 = legFrame(4,3);
t882 = sin(t897);
t886 = cos(t897);
t767 = (-t830 * t886 + t834 * t882) * t859 + t858 * (t830 * t882 + t834 * t886);
t1049 = t767 * t1069;
t893 = pkin(7) + qJ(2,3);
t879 = qJ(3,3) + t893;
t864 = sin(t879);
t867 = cos(t879);
t873 = sin(t893);
t876 = cos(t893);
t819 = t864 * t876 - t867 * t873;
t1068 = 0.1e1 / t819;
t914 = koppelP(3,2);
t918 = koppelP(3,1);
t843 = t890 * t918 + t891 * t914;
t831 = t843 * t909 - t911;
t847 = -t890 * t914 + t891 * t918;
t835 = t847 * t909 + t910;
t898 = legFrame(3,3);
t883 = sin(t898);
t887 = cos(t898);
t771 = (-t831 * t887 + t835 * t883) * t867 + t864 * (t831 * t883 + t835 * t887);
t1047 = t771 * t1068;
t894 = pkin(7) + qJ(2,2);
t880 = qJ(3,2) + t894;
t865 = sin(t880);
t868 = cos(t880);
t874 = sin(t894);
t877 = cos(t894);
t820 = t865 * t877 - t868 * t874;
t1067 = 0.1e1 / t820;
t915 = koppelP(2,2);
t919 = koppelP(2,1);
t844 = t890 * t919 + t891 * t915;
t832 = t844 * t909 - t911;
t848 = -t890 * t915 + t891 * t919;
t836 = t848 * t909 + t910;
t899 = legFrame(2,3);
t884 = sin(t899);
t888 = cos(t899);
t772 = (-t832 * t888 + t836 * t884) * t868 + t865 * (t832 * t884 + t836 * t888);
t1045 = t772 * t1067;
t895 = pkin(7) + qJ(2,1);
t881 = qJ(3,1) + t895;
t866 = sin(t881);
t869 = cos(t881);
t875 = sin(t895);
t878 = cos(t895);
t821 = t866 * t878 - t869 * t875;
t1066 = 0.1e1 / t821;
t916 = koppelP(1,2);
t920 = koppelP(1,1);
t845 = t890 * t920 + t891 * t916;
t833 = t845 * t909 - t911;
t849 = -t890 * t916 + t891 * t920;
t837 = t849 * t909 + t910;
t900 = legFrame(1,3);
t885 = sin(t900);
t889 = cos(t900);
t773 = (-t833 * t889 + t837 * t885) * t869 + t866 * (t833 * t885 + t837 * t889);
t1043 = t773 * t1066;
t1073 = -t866 * t885 + t889 * t869;
t1072 = -t865 * t884 + t888 * t868;
t1071 = -t864 * t883 + t887 * t867;
t1070 = -t858 * t882 + t886 * t859;
t923 = 0.1e1 / pkin(2);
t1037 = t1069 * t923;
t922 = 0.1e1 / pkin(3);
t1014 = t922 * t923;
t783 = -pkin(2) * (t870 * t882 - t871 * t886) + t1070 * pkin(3);
t774 = t783 * t830 * t1069 * t1014;
t823 = t858 * t886 + t859 * t882;
t782 = pkin(2) * (t870 * t886 + t871 * t882) + t823 * pkin(3);
t984 = t834 * t782 * t922;
t748 = t774 + (t767 - t984) * t1037;
t956 = t858 * t870 + t859 * t871;
t924 = 0.1e1 / pkin(2) ^ 2;
t991 = t924 * t1049;
t976 = (pkin(3) * t748 + t956 * t1049) * t991;
t1034 = t1068 * t923;
t788 = -pkin(2) * (t873 * t883 - t876 * t887) + t1071 * pkin(3);
t775 = t788 * t831 * t1068 * t1014;
t827 = t864 * t887 + t867 * t883;
t785 = pkin(2) * (t873 * t887 + t876 * t883) + t827 * pkin(3);
t983 = t835 * t785 * t922;
t753 = t775 + (t771 - t983) * t1034;
t955 = t864 * t873 + t867 * t876;
t987 = t924 * t1047;
t975 = (pkin(3) * t753 + t955 * t1047) * t987;
t1033 = t1067 * t923;
t790 = -pkin(2) * (t874 * t884 - t877 * t888) + t1072 * pkin(3);
t776 = t790 * t832 * t1067 * t1014;
t828 = t865 * t888 + t868 * t884;
t786 = pkin(2) * (t874 * t888 + t877 * t884) + t828 * pkin(3);
t982 = t836 * t786 * t922;
t755 = t776 + (t772 - t982) * t1033;
t954 = t865 * t874 + t868 * t877;
t986 = t924 * t1045;
t974 = (pkin(3) * t755 + t954 * t1045) * t986;
t1032 = t1066 * t923;
t792 = -pkin(2) * (t875 * t885 - t878 * t889) + t1073 * pkin(3);
t777 = t792 * t833 * t1066 * t1014;
t829 = t866 * t889 + t869 * t885;
t787 = pkin(2) * (t875 * t889 + t878 * t885) + t829 * pkin(3);
t981 = t837 * t787 * t922;
t757 = t777 + (t773 - t981) * t1032;
t953 = t866 * t875 + t869 * t878;
t985 = t924 * t1043;
t973 = (pkin(3) * t757 + t953 * t1043) * t985;
t992 = t767 ^ 2 / t806 ^ 2 * t1069;
t990 = t771 ^ 2 / t819 ^ 2 * t1068;
t989 = t772 ^ 2 / t820 ^ 2 * t1067;
t988 = t773 ^ 2 / t821 ^ 2 * t1066;
t758 = -t984 * t1037 + t774;
t1000 = t748 * t758 * t923;
t972 = t1069 * t1000;
t742 = (t956 * pkin(2) + pkin(3)) * t972;
t1013 = 0.2e1 * pkin(3);
t896 = t909 ^ 2;
t1015 = t896 * t923;
t921 = pkin(3) ^ 2;
t936 = (-(-t748 * t921 + (-t1049 - t956 * (t774 / 0.2e1 + (t767 - t984 / 0.2e1) * t1037) * t1013) * pkin(2)) * t991 + (-t782 * t842 - t783 * t846) * t1015) * t922;
t944 = (-t1070 * t846 - t823 * t842) * t1015;
t940 = t1069 * t944;
t980 = pkin(3) * t1000;
t723 = -t742 + t940 - (t936 - t976 - t980) * t1069;
t1065 = t723 * t1069;
t759 = -t983 * t1034 + t775;
t999 = t753 * t759 * t923;
t971 = t1068 * t999;
t743 = (t955 * pkin(2) + pkin(3)) * t971;
t935 = (-(-t753 * t921 + (-t1047 - t955 * (t775 / 0.2e1 + (t771 - t983 / 0.2e1) * t1034) * t1013) * pkin(2)) * t987 + (-t785 * t843 - t788 * t847) * t1015) * t922;
t943 = (-t1071 * t847 - t827 * t843) * t1015;
t939 = t1068 * t943;
t979 = pkin(3) * t999;
t727 = -t743 + t939 - (t935 - t975 - t979) * t1068;
t1064 = t727 * t1068;
t760 = -t982 * t1033 + t776;
t998 = t755 * t760 * t923;
t970 = t1067 * t998;
t744 = (t954 * pkin(2) + pkin(3)) * t970;
t934 = (-(-t755 * t921 + (-t1045 - t954 * (t776 / 0.2e1 + (t772 - t982 / 0.2e1) * t1033) * t1013) * pkin(2)) * t986 + (-t786 * t844 - t790 * t848) * t1015) * t922;
t942 = (-t1072 * t848 - t828 * t844) * t1015;
t938 = t1067 * t942;
t978 = pkin(3) * t998;
t728 = -t744 + t938 - (t934 - t974 - t978) * t1067;
t1063 = t728 * t1067;
t761 = -t981 * t1032 + t777;
t997 = t757 * t761 * t923;
t969 = t1066 * t997;
t745 = (t953 * pkin(2) + pkin(3)) * t969;
t933 = (-(-t757 * t921 + (-t1043 - t953 * (t777 / 0.2e1 + (t773 - t981 / 0.2e1) * t1032) * t1013) * pkin(2)) * t985 + (-t787 * t845 - t792 * t849) * t1015) * t922;
t941 = (-t1073 * t849 - t829 * t845) * t1015;
t937 = t1066 * t941;
t977 = pkin(3) * t997;
t729 = -t745 + t937 - (t933 - t973 - t977) * t1066;
t1062 = t729 * t1066;
t730 = pkin(3) * t972 + (t944 + t976) * t1069;
t1061 = t730 * t1069;
t731 = pkin(3) * t971 + (t943 + t975) * t1068;
t1060 = t731 * t1068;
t732 = pkin(3) * t970 + (t942 + t974) * t1067;
t1059 = t732 * t1067;
t733 = pkin(3) * t969 + (t941 + t973) * t1066;
t1058 = t733 * t1066;
t1057 = (t774 + (0.2e1 * t767 - t984) * t1037) * t758;
t1056 = (t775 + (0.2e1 * t771 - t983) * t1034) * t759;
t1055 = (t776 + (0.2e1 * t772 - t982) * t1033) * t760;
t1054 = (t777 + (0.2e1 * t773 - t981) * t1032) * t761;
t794 = t842 * t886 - t846 * t882;
t798 = t842 * t882 + t846 * t886;
t960 = t794 * t859 - t798 * t858;
t762 = pkin(2) * (t794 * t871 - t798 * t870) + t960 * pkin(3);
t1053 = t762 * t1069;
t795 = t843 * t887 - t847 * t883;
t799 = t843 * t883 + t847 * t887;
t959 = t795 * t867 - t799 * t864;
t763 = pkin(2) * (t795 * t876 - t799 * t873) + t959 * pkin(3);
t1052 = t763 * t1068;
t796 = t844 * t888 - t848 * t884;
t800 = t844 * t884 + t848 * t888;
t958 = t796 * t868 - t800 * t865;
t764 = pkin(2) * (t796 * t877 - t800 * t874) + t958 * pkin(3);
t1051 = t764 * t1067;
t797 = t845 * t889 - t849 * t885;
t801 = t845 * t885 + t849 * t889;
t957 = t797 * t869 - t801 * t866;
t765 = pkin(2) * (t797 * t878 - t801 * t875) + t957 * pkin(3);
t1050 = t765 * t1066;
t1041 = t960 * t1069;
t1040 = t959 * t1068;
t1039 = t958 * t1067;
t1038 = t957 * t1066;
t1036 = t1069 * t1070;
t1035 = t1069 * t823;
t1031 = t1068 * t1071;
t1030 = t1068 * t827;
t1029 = t1067 * t1072;
t1028 = t1067 * t828;
t1027 = t1066 * t1073;
t1026 = t1066 * t829;
t1017 = t890 * t896;
t1016 = t891 * t896;
t1012 = t730 * t1053;
t901 = sin(qJ(3,4));
t1011 = t901 * t1061;
t902 = cos(qJ(3,4));
t1010 = t902 * t1061;
t1009 = t731 * t1052;
t903 = sin(qJ(3,3));
t1008 = t903 * t1060;
t906 = cos(qJ(3,3));
t1007 = t906 * t1060;
t1006 = t732 * t1051;
t904 = sin(qJ(3,2));
t1005 = t904 * t1059;
t907 = cos(qJ(3,2));
t1004 = t907 * t1059;
t1003 = t733 * t1050;
t905 = sin(qJ(3,1));
t1002 = t905 * t1058;
t908 = cos(qJ(3,1));
t1001 = t908 * t1058;
t996 = t762 * t992;
t995 = t763 * t990;
t994 = t764 * t989;
t993 = t765 * t988;
t968 = t901 * t992;
t967 = t902 * t992;
t966 = t903 * t990;
t965 = t906 * t990;
t964 = t904 * t989;
t963 = t907 * t989;
t962 = t905 * t988;
t961 = t908 * t988;
t722 = -t742 + 0.2e1 * t940 - (t936 - 0.2e1 * t976 - 0.2e1 * t980) * t1069;
t952 = t902 * t1057 + t722 * t901;
t951 = t901 * t1057 - t722 * t902;
t724 = -t743 + 0.2e1 * t939 - (t935 - 0.2e1 * t975 - 0.2e1 * t979) * t1068;
t950 = t906 * t1056 + t724 * t903;
t949 = t903 * t1056 - t724 * t906;
t725 = -t744 + 0.2e1 * t938 - (t934 - 0.2e1 * t974 - 0.2e1 * t978) * t1067;
t948 = t907 * t1055 + t725 * t904;
t947 = t904 * t1055 - t725 * t907;
t726 = -t745 + 0.2e1 * t937 - (t933 - 0.2e1 * t973 - 0.2e1 * t977) * t1066;
t946 = t908 * t1054 + t726 * t905;
t945 = t905 * t1054 - t726 * t908;
t1 = [0, (t733 * t1027 + t732 * t1029 + t731 * t1031 + t730 * t1036) * t923, 0, 0, (t723 * t1036 + t727 * t1031 + t728 * t1029 + t729 * t1027 + (-t792 * t1062 - t790 * t1063 - t788 * t1064 - t783 * t1065) * t922) * t923, -t945 * t1027 - t947 * t1029 - t949 * t1031 - t951 * t1036 + (-t783 * t1010 - t788 * t1007 - t790 * t1004 - t792 * t1001 + (-t783 * t968 - t788 * t966 - t790 * t964 - t792 * t962) * t924) * t922, -t946 * t1027 - t948 * t1029 - t950 * t1031 - t952 * t1036 + (t783 * t1011 + t788 * t1008 + t790 * t1005 + t792 * t1002 + (-t783 * t967 - t788 * t965 - t790 * t963 - t792 * t961) * t924) * t922, 0, -t1016, t1017, 0; 0, (t733 * t1026 + t732 * t1028 + t731 * t1030 + t730 * t1035) * t923, 0, 0, (t723 * t1035 + t727 * t1030 + t728 * t1028 + t729 * t1026 + (-t787 * t1062 - t786 * t1063 - t785 * t1064 - t782 * t1065) * t922) * t923, -t945 * t1026 - t947 * t1028 - t949 * t1030 - t951 * t1035 + (-t782 * t1010 - t785 * t1007 - t786 * t1004 - t787 * t1001 + (-t782 * t968 - t785 * t966 - t786 * t964 - t787 * t962) * t924) * t922, -t946 * t1026 - t948 * t1028 - t950 * t1030 - t952 * t1035 + (t782 * t1011 + t785 * t1008 + t786 * t1005 + t787 * t1002 + (-t782 * t967 - t785 * t965 - t786 * t963 - t787 * t961) * t924) * t922, 0, -t1017, -t1016, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, (-t733 * t1038 - t732 * t1039 - t731 * t1040 - t730 * t1041) * t923, 0, 0, (-t723 * t1041 - t727 * t1040 - t728 * t1039 - t729 * t1038 + (t729 * t1050 + t728 * t1051 + t727 * t1052 + t723 * t1053) * t922) * t923, t945 * t1038 + t947 * t1039 + t949 * t1040 + t951 * t1041 + (t902 * t1012 + t906 * t1009 + t907 * t1006 + t908 * t1003 + (t901 * t996 + t903 * t995 + t904 * t994 + t905 * t993) * t924) * t922, t946 * t1038 + t948 * t1039 + t950 * t1040 + t952 * t1041 + (-t901 * t1012 - t903 * t1009 - t904 * t1006 - t905 * t1003 + (t902 * t996 + t906 * t995 + t907 * t994 + t908 * t993) * t924) * t922, 0, 0, 0, 0;];
tau_reg  = t1;
