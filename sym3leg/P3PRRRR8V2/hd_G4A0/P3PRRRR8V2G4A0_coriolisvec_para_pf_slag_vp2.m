% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V2G4A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:14:57
% EndTime: 2020-08-06 18:15:06
% DurationCPUTime: 8.27s
% Computational Cost: add. (52884->364), mult. (108027->709), div. (4716->10), fcn. (126546->34), ass. (0->292)
t1027 = cos(qJ(3,1));
t1006 = cos(pkin(4));
t1021 = sin(qJ(3,1));
t1082 = t1006 * t1021;
t1004 = sin(pkin(4));
t1022 = sin(qJ(2,1));
t1093 = t1004 * t1022;
t1002 = t1027 ^ 2;
t1143 = pkin(3) * t1002;
t1028 = cos(qJ(2,1));
t1034 = pkin(7) + pkin(6);
t968 = t1028 * t1034;
t947 = pkin(2) * t1022 - t968;
t917 = pkin(3) * t1082 + t947 * t1004;
t887 = 0.1e1 / (pkin(2) * t1082 + t917 * t1027 + t1093 * t1143);
t1025 = cos(qJ(3,2));
t1019 = sin(qJ(3,2));
t1084 = t1006 * t1019;
t1020 = sin(qJ(2,2));
t1095 = t1004 * t1020;
t1001 = t1025 ^ 2;
t1144 = pkin(3) * t1001;
t1026 = cos(qJ(2,2));
t967 = t1026 * t1034;
t946 = pkin(2) * t1020 - t967;
t916 = pkin(3) * t1084 + t946 * t1004;
t886 = 0.1e1 / (pkin(2) * t1084 + t916 * t1025 + t1095 * t1144);
t1023 = cos(qJ(3,3));
t1017 = sin(qJ(3,3));
t1086 = t1006 * t1017;
t1018 = sin(qJ(2,3));
t1097 = t1004 * t1018;
t1000 = t1023 ^ 2;
t1145 = pkin(3) * t1000;
t1024 = cos(qJ(2,3));
t966 = t1024 * t1034;
t945 = pkin(2) * t1018 - t966;
t915 = pkin(3) * t1086 + t945 * t1004;
t885 = 0.1e1 / (pkin(2) * t1086 + t915 * t1023 + t1097 * t1145);
t1003 = sin(pkin(8));
t1005 = cos(pkin(8));
t1007 = legFrame(3,3);
t971 = sin(t1007);
t977 = cos(t1007);
t924 = t1003 * t977 + t1005 * t971;
t1160 = t1004 * t924;
t1008 = legFrame(2,3);
t972 = sin(t1008);
t978 = cos(t1008);
t925 = t1003 * t978 + t1005 * t972;
t1159 = t1004 * t925;
t1009 = legFrame(1,3);
t973 = sin(t1009);
t979 = cos(t1009);
t926 = t1003 * t979 + t1005 * t973;
t1158 = t1004 * t926;
t1058 = t1023 * mrSges(3,1) - mrSges(3,2) * t1017;
t1057 = t1025 * mrSges(3,1) - mrSges(3,2) * t1019;
t1056 = t1027 * mrSges(3,1) - mrSges(3,2) * t1021;
t1013 = Ifges(3,1) - Ifges(3,2);
t1029 = mrSges(3,2) * pkin(2);
t1151 = -2 * Ifges(3,4);
t1154 = Ifges(3,4) + t1002 * t1151 + (-t1013 * t1021 + t1029) * t1027;
t1153 = Ifges(3,4) + t1001 * t1151 + (-t1013 * t1019 + t1029) * t1025;
t1152 = Ifges(3,4) + t1000 * t1151 + (-t1013 * t1017 + t1029) * t1023;
t1030 = mrSges(3,1) * pkin(2);
t983 = pkin(6) * mrSges(3,2) - Ifges(3,6);
t1150 = -t983 / 0.2e1;
t984 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t1149 = t984 / 0.2e1;
t1031 = xDP(3);
t1032 = xDP(2);
t1038 = 0.1e1 / pkin(3);
t1033 = xDP(1);
t1014 = legFrame(3,2);
t991 = cos(t1014);
t1105 = t1033 * t991;
t1075 = t1018 * t1034;
t963 = t1023 * pkin(3) + pkin(2);
t1111 = t1006 * (t963 * t1024 + t1075);
t1010 = legFrame(3,1);
t974 = sin(t1010);
t988 = sin(t1014);
t1130 = t974 * t988;
t921 = -t1003 * t971 + t977 * t1005;
t980 = cos(t1010);
t882 = t921 * t1130 + t924 * t980;
t939 = t1018 * t963 - t966;
t858 = -t882 * t1111 + (t924 * t1130 - t921 * t980) * t939;
t1127 = t980 * t988;
t879 = t921 * t1127 - t974 * t924;
t861 = t879 * t1111 - (t924 * t1127 + t974 * t921) * t939;
t873 = -t921 * t1111 + t924 * t939;
t1092 = t1004 * t1023;
t942 = t963 * t1086;
t903 = 0.1e1 / (t939 * t1092 + t942);
t840 = (t1031 * t861 + t1032 * t858 + t873 * t1105) * t903 * t1038;
t1148 = pkin(3) * t840;
t1015 = legFrame(2,2);
t992 = cos(t1015);
t1104 = t1033 * t992;
t1073 = t1020 * t1034;
t964 = t1025 * pkin(3) + pkin(2);
t1110 = t1006 * (t964 * t1026 + t1073);
t1011 = legFrame(2,1);
t975 = sin(t1011);
t989 = sin(t1015);
t1129 = t975 * t989;
t922 = -t1003 * t972 + t978 * t1005;
t981 = cos(t1011);
t883 = t922 * t1129 + t925 * t981;
t940 = t1020 * t964 - t967;
t859 = -t883 * t1110 + (t925 * t1129 - t922 * t981) * t940;
t1126 = t981 * t989;
t880 = t922 * t1126 - t975 * t925;
t862 = t880 * t1110 - (t925 * t1126 + t975 * t922) * t940;
t874 = -t922 * t1110 + t925 * t940;
t1090 = t1004 * t1025;
t943 = t964 * t1084;
t904 = 0.1e1 / (t940 * t1090 + t943);
t841 = (t1031 * t862 + t1032 * t859 + t874 * t1104) * t904 * t1038;
t1147 = pkin(3) * t841;
t1016 = legFrame(1,2);
t993 = cos(t1016);
t1103 = t1033 * t993;
t1071 = t1022 * t1034;
t965 = t1027 * pkin(3) + pkin(2);
t1109 = t1006 * (t965 * t1028 + t1071);
t1012 = legFrame(1,1);
t976 = sin(t1012);
t990 = sin(t1016);
t1128 = t976 * t990;
t923 = -t1003 * t973 + t979 * t1005;
t982 = cos(t1012);
t884 = t923 * t1128 + t926 * t982;
t941 = t1022 * t965 - t968;
t860 = -t884 * t1109 + (t926 * t1128 - t923 * t982) * t941;
t1125 = t982 * t990;
t881 = t923 * t1125 - t976 * t926;
t863 = t881 * t1109 - (t926 * t1125 + t976 * t923) * t941;
t875 = -t923 * t1109 + t926 * t941;
t1088 = t1004 * t1027;
t944 = t965 * t1082;
t905 = 0.1e1 / (t941 * t1088 + t944);
t842 = (t1031 * t863 + t1032 * t860 + t875 * t1103) * t905 * t1038;
t1146 = pkin(3) * t842;
t1142 = t1017 * pkin(2);
t1141 = t1019 * pkin(2);
t1140 = t1021 * pkin(2);
t1136 = t1017 * mrSges(3,1);
t1135 = t1019 * mrSges(3,1);
t1134 = t1021 * mrSges(3,1);
t1037 = pkin(3) ^ 2;
t1067 = t1017 * t1148;
t1124 = 0.2e1 * pkin(2) * pkin(3);
t1085 = t1006 * t1018;
t927 = t1003 * t1085 - t1005 * t1024;
t930 = t1003 * t1024 + t1005 * t1085;
t1046 = t927 * t977 + t971 * t930;
t888 = -t971 * t927 + t930 * t977;
t852 = (-t974 * t1046 + t888 * t1127) * t1017 + t879 * t1092;
t855 = (-t1046 * t980 - t888 * t1130) * t1017 - t882 * t1092;
t870 = t921 * t1092 + t1017 * (t1024 * t924 + t921 * t1085);
t1076 = t1018 * t1023;
t918 = pkin(3) * t1076 + t945;
t900 = 0.1e1 / (t918 * t1092 + t942);
t834 = -t870 * t885 * t1105 + (t1031 * t852 + t1032 * t855) * t900;
t1039 = pkin(2) ^ 2;
t969 = t1034 ^ 2 + t1039;
t1133 = (-t1034 * t1067 + (t1000 * t1037 + t1023 * t1124 + t969) * t834) * t834;
t1066 = t1019 * t1147;
t1083 = t1006 * t1020;
t928 = t1003 * t1083 - t1005 * t1026;
t931 = t1003 * t1026 + t1005 * t1083;
t1045 = t928 * t978 + t972 * t931;
t889 = -t972 * t928 + t931 * t978;
t853 = (-t975 * t1045 + t889 * t1126) * t1019 + t880 * t1090;
t856 = (-t1045 * t981 - t889 * t1129) * t1019 - t883 * t1090;
t871 = t922 * t1090 + t1019 * (t1026 * t925 + t922 * t1083);
t1074 = t1020 * t1025;
t919 = pkin(3) * t1074 + t946;
t901 = 0.1e1 / (t919 * t1090 + t943);
t835 = -t871 * t886 * t1104 + (t1031 * t853 + t1032 * t856) * t901;
t1132 = (-t1034 * t1066 + (t1001 * t1037 + t1025 * t1124 + t969) * t835) * t835;
t1065 = t1021 * t1146;
t1081 = t1006 * t1022;
t929 = t1003 * t1081 - t1005 * t1028;
t932 = t1003 * t1028 + t1005 * t1081;
t1044 = t929 * t979 + t973 * t932;
t890 = -t973 * t929 + t932 * t979;
t854 = (-t976 * t1044 + t890 * t1125) * t1021 + t881 * t1088;
t857 = (-t1044 * t982 - t890 * t1128) * t1021 - t884 * t1088;
t872 = t923 * t1088 + t1021 * (t1028 * t926 + t923 * t1081);
t1072 = t1022 * t1027;
t920 = pkin(3) * t1072 + t947;
t902 = 0.1e1 / (t920 * t1088 + t944);
t836 = -t872 * t887 * t1103 + (t1031 * t854 + t1032 * t857) * t902;
t1131 = (-t1034 * t1065 + (t1002 * t1037 + t1027 * t1124 + t969) * t836) * t836;
t994 = m(3) * pkin(2) + mrSges(2,1);
t1108 = t1024 * t834;
t970 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t1117 = t1004 * ((t1058 + t994) * t1024 + t970 * t1018);
t1091 = t1004 * t1024;
t1098 = t1004 * t1017;
t1102 = t834 * t1034;
t1064 = t1017 * t1102;
t825 = t1064 - t1148;
t807 = (((t1006 * t840 + t834 * t1091) * t1145 + ((-t1067 + t1102) * t1018 + pkin(2) * t1108) * t1092 + t825 * t1006) * t834 + (t840 * t1091 + (t1000 * t1006 - t1076 * t1098 - t1006) * t834) * t1148) * t885;
t1080 = t1006 * t1038;
t810 = t885 * t1080 * t1133 + (-t1006 * t1064 + (-t918 * t1098 + (pkin(2) * t1023 + t1145) * t1006) * t840) * t900 * t840;
t816 = (-t1023 * t1133 - (pkin(2) * t840 - t825 * t1023) * t1148) * t885;
t831 = t834 ^ 2;
t837 = t840 ^ 2;
t954 = t1023 * mrSges(3,2) + t1136;
t906 = -t1058 * t1006 + t954 * t1097;
t998 = -m(1) - m(2) - m(3);
t1123 = -t807 * t1117 + t906 * t810 + t998 * t816 + ((-t831 * t994 - t1058 * (t831 + t837)) * t1018 + (t834 * t970 - 0.2e1 * t840 * t954) * t1108) * t1004 - t837 * t1006 * t954;
t1107 = t1026 * t835;
t1116 = t1004 * ((t1057 + t994) * t1026 + t970 * t1020);
t1089 = t1004 * t1026;
t1096 = t1004 * t1019;
t1101 = t835 * t1034;
t1063 = t1019 * t1101;
t826 = t1063 - t1147;
t808 = (((t1006 * t841 + t835 * t1089) * t1144 + ((-t1066 + t1101) * t1020 + pkin(2) * t1107) * t1090 + t826 * t1006) * t835 + (t841 * t1089 + (t1001 * t1006 - t1074 * t1096 - t1006) * t835) * t1147) * t886;
t811 = t886 * t1080 * t1132 + (-t1006 * t1063 + (-t919 * t1096 + (pkin(2) * t1025 + t1144) * t1006) * t841) * t901 * t841;
t817 = (-t1025 * t1132 - (pkin(2) * t841 - t826 * t1025) * t1147) * t886;
t832 = t835 ^ 2;
t838 = t841 ^ 2;
t955 = t1025 * mrSges(3,2) + t1135;
t907 = -t1057 * t1006 + t955 * t1095;
t1122 = -t808 * t1116 + t907 * t811 + t998 * t817 + ((-t832 * t994 - t1057 * (t832 + t838)) * t1020 + (t835 * t970 - 0.2e1 * t841 * t955) * t1107) * t1004 - t838 * t1006 * t955;
t1106 = t1028 * t836;
t1115 = t1004 * ((t1056 + t994) * t1028 + t970 * t1022);
t1087 = t1004 * t1028;
t1094 = t1004 * t1021;
t1100 = t836 * t1034;
t1062 = t1021 * t1100;
t827 = t1062 - t1146;
t809 = (((t1006 * t842 + t836 * t1087) * t1143 + ((-t1065 + t1100) * t1022 + pkin(2) * t1106) * t1088 + t827 * t1006) * t836 + (t842 * t1087 + (t1002 * t1006 - t1072 * t1094 - t1006) * t836) * t1146) * t887;
t812 = t887 * t1080 * t1131 + (-t1006 * t1062 + (-t920 * t1094 + (pkin(2) * t1027 + t1143) * t1006) * t842) * t902 * t842;
t818 = (-t1027 * t1131 - (pkin(2) * t842 - t827 * t1027) * t1146) * t887;
t833 = t836 ^ 2;
t839 = t842 ^ 2;
t956 = t1027 * mrSges(3,2) + t1134;
t908 = -t1056 * t1006 + t956 * t1093;
t1121 = -t809 * t1115 + t908 * t812 + t998 * t818 + ((-t833 * t994 - t1056 * (t833 + t839)) * t1022 + (t836 * t970 - 0.2e1 * t842 * t956) * t1106) * t1004 - t839 * t1006 * t956;
t1114 = t1004 * t921;
t1113 = t1004 * t922;
t1112 = t1004 * t923;
t1099 = 0.2e1 * t1029;
t1061 = t1123 * t885;
t1060 = t1122 * t886;
t1059 = t1121 * t887;
t1040 = -0.2e1 * pkin(6) * mrSges(3,3) + (-pkin(6) ^ 2 - t1039) * m(3) - Ifges(3,1) - Ifges(2,3);
t933 = t1017 * t984 + t983 * t1023;
t1055 = 0.2e1 * ((t1150 * t1017 + t1149 * t1023) * t840 + (t1030 * t1017 + t1152) * t834) * t840 + t816 * t1117 - (t1013 * t1000 - 0.2e1 * (Ifges(3,4) * t1017 + t1030) * t1023 + t1017 * t1099 + t1040) * t807 - t933 * t810;
t934 = t1019 * t984 + t983 * t1025;
t1054 = 0.2e1 * ((t1150 * t1019 + t1149 * t1025) * t841 + (t1030 * t1019 + t1153) * t835) * t841 + t817 * t1116 - (t1013 * t1001 - 0.2e1 * (Ifges(3,4) * t1019 + t1030) * t1025 + t1019 * t1099 + t1040) * t808 - t934 * t811;
t935 = t1021 * t984 + t983 * t1027;
t1053 = 0.2e1 * ((t1150 * t1021 + t1149 * t1027) * t842 + (t1030 * t1021 + t1154) * t836) * t842 + t818 * t1115 - (t1013 * t1002 - 0.2e1 * (Ifges(3,4) * t1021 + t1030) * t1027 + t1021 * t1099 + t1040) * t809 - t935 * t812;
t1052 = t1055 * t900;
t1051 = t1054 * t901;
t1050 = t1053 * t902;
t1049 = (-Ifges(3,3) * t810 + t933 * t807 + t906 * t816 + t831 * (pkin(2) * t1136 + t1152)) * t903;
t1048 = (-Ifges(3,3) * t811 + t934 * t808 + t907 * t817 + t832 * (pkin(2) * t1135 + t1153)) * t904;
t1047 = (-Ifges(3,3) * t812 + t935 * t809 + t908 * t818 + t833 * (pkin(2) * t1134 + t1154)) * t905;
t1043 = pkin(3) * t1098 - t1006 * t945;
t1042 = pkin(3) * t1096 - t1006 * t946;
t1041 = pkin(3) * t1094 - t1006 * t947;
t950 = pkin(2) * t1028 + t1071;
t949 = pkin(2) * t1026 + t1073;
t948 = pkin(2) * t1024 + t1075;
t899 = t1006 * t993 - t990 * t1158;
t898 = t1006 * t992 - t989 * t1159;
t897 = t1006 * t991 - t988 * t1160;
t896 = t950 * t1003 - t1041 * t1005;
t895 = t949 * t1003 - t1042 * t1005;
t894 = t948 * t1003 - t1043 * t1005;
t893 = -t1041 * t1003 - t950 * t1005;
t892 = -t1042 * t1003 - t949 * t1005;
t891 = -t1043 * t1003 - t948 * t1005;
t878 = t1044 * t990 + t993 * t1093;
t877 = t1045 * t989 + t992 * t1095;
t876 = t1046 * t988 + t991 * t1097;
t869 = -t973 * t893 + t896 * t979;
t868 = -t972 * t892 + t895 * t978;
t867 = -t971 * t891 + t894 * t977;
t866 = t917 * t993 + (t893 * t979 + t896 * t973) * t990;
t865 = t916 * t992 + (t892 * t978 + t895 * t972) * t989;
t864 = t915 * t991 + (t891 * t977 + t894 * t971) * t988;
t1 = [(t1053 * t993 * t872 + t1121 * (-((-t1028 * t923 + t926 * t1081) * t993 - t990 * t1093) * t1143 + ((t1041 * t926 + t950 * t923) * t993 + t917 * t990) * t1027 + (t1006 * t990 + t993 * t1158) * t1140)) * t887 + (t1054 * t992 * t871 + t1122 * (-((-t1026 * t922 + t925 * t1083) * t992 - t989 * t1095) * t1144 + ((t1042 * t925 + t949 * t922) * t992 + t916 * t989) * t1025 + (t1006 * t989 + t992 * t1159) * t1141)) * t886 + (t1055 * t991 * t870 + t1123 * (-((-t1024 * t921 + t924 * t1085) * t991 - t988 * t1097) * t1145 + ((t1043 * t924 + t948 * t921) * t991 + t915 * t988) * t1023 + (t1006 * t988 + t991 * t1160) * t1142)) * t885 + (t993 * t875 * t1047 + t992 * t874 * t1048 + t991 * t873 * t1049) * t1038; -t857 * t1050 - t856 * t1051 - t855 * t1052 + (-(t878 * t976 - t890 * t982) * t1143 + (-t866 * t976 + t869 * t982) * t1027 - (t982 * t1112 + t899 * t976) * t1140) * t1059 + (-(t877 * t975 - t889 * t981) * t1144 + (-t865 * t975 + t868 * t981) * t1025 - (t981 * t1113 + t898 * t975) * t1141) * t1060 + (-(t876 * t974 - t888 * t980) * t1145 + (-t864 * t974 + t867 * t980) * t1023 - (t980 * t1114 + t897 * t974) * t1142) * t1061 + (t860 * t1047 + t859 * t1048 + t858 * t1049) * t1038; -t854 * t1050 - t853 * t1051 - t852 * t1052 + ((t878 * t982 + t890 * t976) * t1143 + (t866 * t982 + t976 * t869) * t1027 + (-t1112 * t976 + t899 * t982) * t1140) * t1059 + ((t877 * t981 + t889 * t975) * t1144 + (t865 * t981 + t975 * t868) * t1025 + (-t1113 * t975 + t898 * t981) * t1141) * t1060 + ((t876 * t980 + t888 * t974) * t1145 + (t864 * t980 + t974 * t867) * t1023 + (-t1114 * t974 + t897 * t980) * t1142) * t1061 + (t1047 * t863 + t1048 * t862 + t1049 * t861) * t1038;];
taucX  = t1;
