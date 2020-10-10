% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR6V1G2A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:36:19
% EndTime: 2020-08-06 18:36:27
% DurationCPUTime: 7.47s
% Computational Cost: add. (43605->450), mult. (48396->647), div. (7563->16), fcn. (33165->102), ass. (0->326)
t963 = sin(pkin(7));
t1155 = t963 * pkin(1);
t1190 = pkin(5) + t1155;
t1199 = t1190 * mrSges(3,2);
t964 = cos(pkin(7));
t902 = t964 * pkin(1);
t874 = t902 + pkin(2);
t1198 = mrSges(3,1) * t874;
t983 = cos(qJ(1,1));
t990 = -pkin(6) - pkin(5);
t1117 = t983 * t990;
t1118 = t983 * t963;
t1114 = t990 * t963;
t863 = -pkin(1) + t1114;
t977 = sin(qJ(1,1));
t846 = t863 * t977;
t1197 = pkin(2) * t1118 + (t977 * pkin(2) + t1117) * t964 - t846;
t981 = cos(qJ(1,2));
t1119 = t981 * t990;
t1120 = t981 * t963;
t975 = sin(qJ(1,2));
t845 = t863 * t975;
t1196 = pkin(2) * t1120 + (t975 * pkin(2) + t1119) * t964 - t845;
t979 = cos(qJ(1,3));
t1121 = t979 * t990;
t1122 = t979 * t963;
t973 = sin(qJ(1,3));
t844 = t863 * t973;
t1195 = pkin(2) * t1122 + (t973 * pkin(2) + t1121) * t964 - t844;
t1165 = -2 * Ifges(3,4);
t1191 = mrSges(3,2) * t874;
t982 = cos(qJ(3,1));
t961 = t982 ^ 2;
t968 = Ifges(3,1) - Ifges(3,2);
t976 = sin(qJ(3,1));
t1194 = t961 * t1165 + Ifges(3,4) + (-t968 * t976 + t1191) * t982;
t980 = cos(qJ(3,2));
t960 = t980 ^ 2;
t974 = sin(qJ(3,2));
t1193 = t960 * t1165 + Ifges(3,4) + (-t968 * t974 + t1191) * t980;
t978 = cos(qJ(3,3));
t959 = t978 ^ 2;
t972 = sin(qJ(3,3));
t1192 = t959 * t1165 + Ifges(3,4) + (-t968 * t972 + t1191) * t978;
t1012 = 0.1e1 / pkin(3);
t1154 = t974 * pkin(2);
t1092 = 0.2e1 * t1154;
t1005 = 0.2e1 * qJ(3,2);
t942 = sin(t1005);
t1156 = pkin(3) * t942;
t1162 = 0.2e1 * t990;
t1167 = -0.2e1 * pkin(2);
t1168 = -0.2e1 * pkin(1);
t923 = pkin(7) + qJ(3,2);
t889 = sin(t923);
t927 = -pkin(7) + qJ(3,2);
t891 = sin(t927);
t935 = qJ(1,2) + pkin(7);
t893 = sin(t935);
t1141 = (t893 * t1162 + cos(t935) * t1167 + t981 * t1168 + (-cos(qJ(1,2) - t927) - cos(qJ(1,2) + t923)) * pkin(3)) / (t1092 + t1156 + (t889 + t891) * pkin(1));
t987 = xDP(3);
t1066 = t987 * t1141;
t915 = t980 * pkin(3);
t883 = t915 + pkin(2);
t1136 = t883 * t975;
t848 = 0.1e1 / (t902 + t883);
t957 = 0.1e1 / t974;
t1179 = t848 * t957;
t1068 = ((t1119 + t1136) * t964 - t845 + t883 * t1120) * t1179;
t970 = legFrame(2,2);
t912 = cos(t970);
t989 = xDP(1);
t1132 = t912 * t989;
t909 = sin(t970);
t988 = xDP(2);
t1133 = t909 * t988;
t1189 = (t1066 / 0.6e1 - (t1132 / 0.6e1 - t1133 / 0.6e1) * t1068) * t1012;
t1011 = pkin(3) ^ 2;
t1188 = (t1066 / 0.4e1 - (t1132 / 0.4e1 - t1133 / 0.4e1) * t1068) * t1011 * t1012;
t1187 = (t1066 / 0.3e1 - (t1132 / 0.3e1 - t1133 / 0.3e1) * t1068) * t1012;
t1186 = (t1066 / 0.2e1 - (t1132 / 0.2e1 - t1133 / 0.2e1) * t1068) * t1012;
t1185 = 0.2e1 * pkin(2);
t1184 = 0.4e1 * pkin(2);
t1116 = t988 / 0.2e1;
t1115 = t989 / 0.2e1;
t1013 = pkin(2) ^ 2;
t1014 = pkin(1) ^ 2;
t992 = m(2) + m(3);
t1017 = -(pkin(5) ^ 2 + t1013) * m(3) - t992 * t1014 - Ifges(3,1) - Ifges(1,3) - Ifges(2,3) - 0.2e1 * (m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3)) * t1155 - 0.2e1 * mrSges(3,3) * pkin(5);
t1032 = Ifges(3,6) / 0.2e1 - t1199 / 0.2e1;
t1084 = mrSges(3,1) * t1155;
t904 = mrSges(3,1) * pkin(5) - Ifges(3,5);
t1033 = t1084 / 0.2e1 + t904 / 0.2e1;
t1070 = -m(3) * pkin(2) - mrSges(2,1);
t1091 = t976 * t1185;
t1094 = 0.2e1 * t902;
t1010 = pkin(3) * t1011;
t962 = t990 ^ 2;
t1104 = t1014 + t962 / 0.2e1;
t1043 = 0.3e1 / 0.8e1 * t1011 + t1013 / 0.2e1 + t1104;
t1082 = pkin(1) * t1114;
t1103 = t1014 + t962;
t999 = 0.2e1 * pkin(7);
t918 = sin(t999);
t919 = cos(t999);
t993 = 0.3e1 * t1013;
t1016 = -0.8e1 * (0.3e1 / 0.4e1 * t1011 + t993 + t1103) * t902 - 0.8e1 * (t1043 - t1082) * t1185 + 0.8e1 * (-pkin(2) * t919 + t918 * t990) * t1014;
t1007 = 0.3e1 * qJ(3,1);
t953 = cos(t1007);
t1026 = pkin(3) * (-t982 * pkin(2) + t874 * t953);
t1105 = t1014 * t918;
t1072 = 0.2e1 * t1105;
t1081 = t990 * t902;
t1027 = t1072 - 0.4e1 * t1081;
t887 = t1014 * t919;
t1174 = -t1011 / 0.2e1 - t887;
t1063 = -0.2e1 * t1013 + t1174;
t1042 = -t1014 + t1063;
t998 = -0.2e1 * t1014;
t1044 = -0.2e1 * t887 - 0.4e1 * t1013 + t998 - t1011;
t1046 = -0.2e1 * t1081 + t1105;
t1099 = 0.2e1 * pkin(3);
t867 = pkin(6) + t1190;
t1059 = t867 * t1099;
t1161 = -0.6e1 * t1011;
t1073 = t874 * t1161;
t1074 = -0.2e1 * t1011 * t867;
t916 = t982 * pkin(3);
t884 = t916 + t1185;
t1083 = t884 * t902;
t1087 = pkin(2) * t902;
t865 = -0.2e1 * t1082;
t997 = 0.2e1 * t1014;
t1095 = -0.4e1 * pkin(3) * (0.6e1 * t1087 + t865 + t997 + t993 + t962 - t1174);
t1098 = 0.4e1 * pkin(3);
t1006 = 0.4e1 * qJ(3,1);
t1106 = t1011 * cos(t1006);
t994 = 0.2e1 * t1013;
t1102 = t994 + t1014;
t840 = 0.4e1 * t1087 + t887 + t1102;
t1008 = 0.2e1 * qJ(3,1);
t954 = cos(t1008);
t1138 = t840 * t954;
t1159 = pkin(2) * t867;
t830 = t1046 + 0.2e1 * t1159;
t1139 = t830 * t954;
t1153 = t1012 / 0.2e1;
t1164 = -0.2e1 * t1011 - 0.2e1 * t840;
t1055 = -t990 + t1155;
t885 = t916 + pkin(2);
t1135 = t885 * t977;
t849 = 0.1e1 / (t902 + t885);
t958 = 0.1e1 / t976;
t1178 = t849 * t958;
t1067 = ((t1117 + t1135) * t964 - t846 + t885 * t1118) * t1178;
t1089 = pkin(7) + qJ(3,1);
t1090 = -pkin(7) + qJ(3,1);
t936 = qJ(1,1) + pkin(7);
t894 = sin(t936);
t945 = sin(t1008);
t1140 = (t894 * t1162 + cos(t936) * t1167 + t983 * t1168 + (-cos(qJ(1,1) - t1090) - cos(qJ(1,1) + t1089)) * pkin(3)) / (t1091 + pkin(3) * t945 + (sin(t1089) + sin(t1090)) * pkin(1));
t971 = legFrame(1,2);
t910 = sin(t971);
t913 = cos(t971);
t787 = (t987 * t1140 + (t910 * t988 - t913 * t989) * t1067) * t1012;
t1147 = t787 * t976;
t1075 = pkin(3) * t1147;
t872 = t971 + t936;
t873 = -t971 + t936;
t836 = -sin(t872) + sin(t873);
t839 = cos(t873) + cos(t872);
t801 = (t839 * t1115 + t836 * t1116 - t894 * t987) * t849;
t775 = t1055 * t801 - t1075;
t943 = sin(t1006);
t944 = sin(t1007);
t765 = ((t953 * t1074 + (t884 * t867 + t1046 - t1139) * t1099) * t787 + (-t1010 * t943 + t1016 * t976 + t944 * t1073 + t945 * t1095) * t801) / (t1138 + t1106 / 0.2e1 - 0.2e1 * t1083 + 0.2e1 * t1026 + t1042) * t801 * t1153 + (t775 * t1184 + (-t1011 * t943 + t945 * t1164 + (-t874 * t944 - t976 * t902) * t1098) * t787 + (-0.2e1 * t1139 + (-t953 + t982) * t1059 + t1027) * t801) / (0.4e1 * t1026 + t1044 - 0.4e1 * t1083 + t1106 + 0.2e1 * t1138) * t787;
t772 = (t775 - t1075) * t849 * t801;
t854 = Ifges(3,6) - t1199;
t855 = t1084 + t904;
t825 = -t854 * t982 + t976 * t855;
t907 = t976 * mrSges(3,2);
t1047 = 0.2e1 * ((t1032 * t976 + t1033 * t982) * t787 + (t1198 * t976 + t1194) * t801) * t787 - (t968 * t961 - 0.2e1 * (Ifges(3,4) * t976 + t1198) * t982 + (t1070 + t907) * t1094 + mrSges(3,2) * t1091 + t1017) * t772 - t825 * t765;
t1182 = -t849 * t1047 / 0.2e1;
t1065 = t893 * t848 * t987;
t870 = t970 + t935;
t871 = -t970 + t935;
t835 = -sin(t870) + sin(t871);
t820 = t835 * t848 * t1116;
t838 = cos(t871) + cos(t870);
t822 = t838 * t848 * t1115;
t800 = t820 + t822 - t1065;
t795 = t1014 * t800;
t1018 = t795 + (0.3e1 / 0.4e1 * t822 + 0.3e1 / 0.4e1 * t820 - 0.3e1 / 0.4e1 * t1065) * t1011 + 0.3e1 * (t1013 + t962 / 0.3e1) * t800;
t1003 = 0.4e1 * qJ(3,2);
t1107 = t1011 * cos(t1003);
t1004 = 0.3e1 * qJ(3,2);
t950 = cos(t1004);
t1109 = t950 - t980;
t786 = (t1066 + (-t1132 + t1133) * t1068) * t1012;
t1157 = pkin(3) * t786;
t1158 = pkin(2) * t990;
t1160 = pkin(2) * pkin(3);
t1163 = -0.2e1 * t964;
t1175 = 0.4e1 * t990;
t1148 = t786 * t974;
t1076 = pkin(3) * t1148;
t777 = t1055 * t800 - t1076;
t778 = t800 - t1186;
t779 = t800 + t1186;
t877 = 0.2e1 * t923;
t879 = 0.2e1 * t927;
t922 = pkin(7) + t1005;
t896 = cos(t922);
t897 = cos(t923);
t926 = -pkin(7) + t1005;
t899 = cos(t926);
t900 = cos(t927);
t931 = t1004 + pkin(7);
t932 = t1004 - pkin(7);
t933 = qJ(3,2) + t999;
t934 = qJ(3,2) - 0.2e1 * pkin(7);
t940 = sin(t1003);
t941 = sin(t1004);
t951 = cos(t1005);
t758 = (t1011 * t786 * t950 * t1162 + t1072 * t1157 + 0.4e1 * (t951 * t1158 + t1159 + (-t902 - t915 / 0.2e1) * t990) * t1157 + (-0.8e1 * (t1011 / 0.4e1 + 0.3e1 / 0.2e1 * t1013 + t1104) * t1156 + pkin(2) * t941 * t1161 - t1010 * t940 - 0.16e2 * t1043 * t1154) * t800 + (t778 * sin(t879) + t779 * sin(t877)) * pkin(3) * t998 + ((cos(t934) - cos(t933)) * t1175 - 0.4e1 * (sin(t934) + sin(t933)) * pkin(2)) * t795 + (-0.4e1 * (t1018 + t1188) * t891 - 0.4e1 * (t1018 - t1188) * t889 + (-t778 * t896 + t779 * t899) * pkin(3) * t1175 - 0.3e1 * ((t800 - t1187) * sin(t932) + (t800 + t1187) * sin(t931)) * t1011 - 0.12e2 * ((t800 - t1189) * sin(t926) + (t800 + t1189) * sin(t922)) * t1160 + 0.8e1 * (-t897 + t900) * t800 * t1158) * pkin(1)) / (t994 * t951 + t1107 / 0.2e1 + (cos(t879) / 0.2e1 + cos(t877) / 0.2e1 + t951 - 0.1e1) * t1014 + t1109 * pkin(2) * t1099 + ((t896 + t899 + t1163) * t1185 + (cos(t932) + cos(t931) - t897 - t900) * pkin(3)) * pkin(1) + t1063) * t800 * t1153 + (t777 * t1184 + (-t1011 * t940 + t942 * t1164 + (-t874 * t941 - t974 * t902) * t1098) * t786 + (-t1109 * t1059 - 0.2e1 * t830 * t951 + t1027) * t800) / (0.2e1 * t840 * t951 + t1107 - 0.4e1 * (t915 + t1185) * t902 + (-t980 * pkin(2) + t874 * t950) * t1098 + t1044) * t786;
t774 = (t777 - t1076) * t848 * t800;
t824 = -t854 * t980 + t974 * t855;
t906 = t974 * mrSges(3,2);
t1051 = 0.2e1 * ((t1032 * t974 + t1033 * t980) * t786 + (t1198 * t974 + t1193) * t800) * t786 - (t968 * t960 - 0.2e1 * (Ifges(3,4) * t974 + t1198) * t980 + (t1070 + t906) * t1094 + mrSges(3,2) * t1092 + t1017) * t774 - t824 * t758;
t1181 = -t848 * t1051 / 0.2e1;
t1093 = t972 * t1185;
t914 = t978 * pkin(3);
t882 = t914 + pkin(2);
t1137 = t882 * t973;
t956 = 0.1e1 / t972;
t1143 = ((t1121 + t1137) * t964 - t844 + t882 * t1122) * t956;
t969 = legFrame(3,2);
t911 = cos(t969);
t1069 = t911 * t1143;
t1034 = t1012 * t989 * t1069;
t1001 = 0.3e1 * qJ(3,3);
t947 = cos(t1001);
t1060 = pkin(3) * (t947 - t978);
t1002 = 0.2e1 * qJ(3,3);
t920 = pkin(7) + t1002;
t924 = -pkin(7) + t1002;
t1064 = cos(t920) + cos(t924) + t1163;
t1000 = 0.4e1 * qJ(3,3);
t1108 = t1011 * cos(t1000);
t930 = qJ(1,3) + pkin(7);
t868 = t969 + t930;
t869 = -t969 + t930;
t834 = -sin(t868) + sin(t869);
t837 = cos(t869) + cos(t868);
t847 = 0.1e1 / (t902 + t882);
t1111 = (t1115 * t837 + t1116 * t834) * t847;
t908 = sin(t969);
t1045 = t847 * t908 * t1143;
t921 = pkin(7) + qJ(3,3);
t888 = sin(t921);
t925 = -pkin(7) + qJ(3,3);
t890 = sin(t925);
t892 = sin(t930);
t939 = sin(t1002);
t1142 = (t892 * t1162 + cos(t930) * t1167 + t979 * t1168 + (-cos(qJ(1,3) - t925) - cos(qJ(1,3) + t921)) * pkin(3)) / (t1093 + pkin(3) * t939 + (t888 + t890) * pkin(1));
t1112 = (t1045 * t988 + t1142 * t987) * t1012;
t1134 = t892 * t987;
t785 = -t847 * t1034 + t1112;
t1149 = t785 * t972;
t1077 = pkin(3) * t1149;
t799 = -t847 * t1134 + t1111;
t776 = t1055 * t799 - t1077;
t1113 = 0.2e1 * t1112;
t780 = (0.2e1 * t1034 - t1134) * t847 + t1111 - t1113;
t781 = (-0.2e1 * t1034 - t1134) * t847 + t1111 + t1113;
t876 = 0.2e1 * t921;
t878 = 0.2e1 * t925;
t881 = t914 + t1185;
t928 = t1001 + pkin(7);
t929 = t1001 - pkin(7);
t937 = sin(t1000);
t938 = sin(t1001);
t948 = cos(t1002);
t761 = ((t947 * t1074 + (-t830 * t948 + t881 * t867 + t1046) * t1099) * t785 + (-t1010 * t937 + t1016 * t972 + t938 * t1073 + t939 * t1095) * t799) / (t840 * t948 + t1108 / 0.2e1 - 0.2e1 * t881 * t902 + (-t978 * pkin(2) + t874 * t947) * t1099 + t1042) * t799 * t1153 + (t776 * t1184 + (((t1034 - t1134) * t847 + t1111 - t1112) * sin(t878) - ((-t1034 - t1134) * t847 + t1111 + t1112) * sin(t876)) * t1014 + (-0.2e1 * (t1011 + t1102) * t939 - 0.4e1 * t938 * t1160 - t1011 * t937) * t785 + (t1072 + (t1184 * t948 + 0.2e1 * t1060) * t990) * t799 + ((sin(t924) * t780 - sin(t920) * t781) * t1185 + t1064 * t799 * t1162 + ((-sin(t928) - t890) * t781 + (sin(t929) + t888) * t780) * pkin(3)) * pkin(1)) / ((t997 + 0.4e1 * t1013) * t948 + t1108 + (cos(t878) + cos(t876)) * t1014 + t1060 * t1184 + (t1064 * t1184 + (cos(t929) + cos(t928) - cos(t925) - cos(t921)) * t1099) * pkin(1) + t1044) * t785;
t773 = (t776 - t1077) * t847 * t799;
t823 = -t854 * t978 + t972 * t855;
t905 = t972 * mrSges(3,2);
t1049 = 0.2e1 * ((t1032 * t972 + t1033 * t978) * t785 + (t1198 * t972 + t1192) * t799) * t785 - (t968 * t959 - 0.2e1 * (Ifges(3,4) * t972 + t1198) * t978 + (t1070 + t905) * t1094 + mrSges(3,2) * t1093 + t1017) * t773 - t823 * t761;
t1180 = -t847 * t1049 / 0.2e1;
t1150 = mrSges(3,1) * t976;
t1144 = t801 * t982;
t784 = t787 ^ 2;
t832 = t1013 + t865 + 0.2e1 * t1087 + t1103;
t771 = ((-t874 - t916) * t784 * pkin(3) + (-(t1011 * t961 + t832) * t801 + (-t874 * t1144 + t867 * t1147) * t1099) * t1144) * t1178;
t856 = -t982 * mrSges(3,1) + t907;
t1048 = t801 ^ 2 * (t874 * t1150 + t1194) - Ifges(3,3) * t765 + t856 * t771 + t825 * t772;
t1177 = t1048 * t1067;
t1151 = mrSges(3,1) * t974;
t1145 = t800 * t980;
t783 = t786 ^ 2;
t770 = ((-t874 - t915) * t783 * pkin(3) + (-(t1011 * t960 + t832) * t800 + (-t874 * t1145 + t867 * t1148) * t1099) * t1145) * t1179;
t858 = -t980 * mrSges(3,1) + t906;
t1052 = t800 ^ 2 * (t874 * t1151 + t1193) - Ifges(3,3) * t758 + t858 * t770 + t824 * t774;
t1176 = t1052 * t1068;
t782 = t785 ^ 2;
t1152 = mrSges(3,1) * t972;
t1146 = t799 * t978;
t1128 = t972 * t908;
t1127 = t972 * t911;
t1126 = t974 * t909;
t1125 = t974 * t912;
t1124 = t976 * t910;
t1123 = t976 * t913;
t1080 = pkin(3) * (t973 * t964 + t1122) * t959;
t1079 = pkin(3) * (t975 * t964 + t1120) * t960;
t1078 = pkin(3) * (t977 * t964 + t1118) * t961;
t769 = ((-t874 - t914) * t782 * pkin(3) + (-(t1011 * t959 + t832) * t799 + (-t874 * t1146 + t867 * t1149) * t1099) * t1146) * t956 * t847;
t857 = -t978 * mrSges(3,1) + t905;
t1050 = t799 ^ 2 * (t874 * t1152 + t1192) - Ifges(3,3) * t761 + t857 * t769 + t823 * t773;
t1037 = (t858 * t758 - t992 * t770 - t783 * (mrSges(3,2) * t980 + t1151)) * t957;
t1036 = (t857 * t761 - t992 * t769 - t782 * (mrSges(3,2) * t978 + t1152)) * t956;
t1035 = (t856 * t765 - t992 * t771 - t784 * (mrSges(3,2) * t982 + t1150)) * t958;
t1025 = t847 * t1036;
t1024 = t848 * t1037;
t1023 = t849 * t1035;
t1 = [(t913 * t1078 + (pkin(3) * t1124 + t1197 * t913) * t982 + t874 * t1124) * t1023 + (t912 * t1079 + (pkin(3) * t1126 + t1196 * t912) * t980 + t874 * t1126) * t1024 + (t911 * t1080 + (pkin(3) * t1128 + t1195 * t911) * t978 + t874 * t1128) * t1025 + t839 * t1182 + t838 * t1181 + t837 * t1180 + (-t1050 * t847 * t1069 - t912 * t1176 - t913 * t1177) * t1012; (-t910 * t1078 + (pkin(3) * t1123 - t1197 * t910) * t982 + t874 * t1123) * t1023 + (-t909 * t1079 + (pkin(3) * t1125 - t1196 * t909) * t980 + t874 * t1125) * t1024 + (-t908 * t1080 + (pkin(3) * t1127 - t1195 * t908) * t978 + t874 * t1127) * t1025 + t836 * t1182 + t835 * t1181 + t834 * t1180 + (t1050 * t1045 + t909 * t1176 + t910 * t1177) * t1012; (t1047 * t894 + t982 * ((t885 * t983 - t977 * t990) * t964 - t863 * t983 - t963 * t1135) * t1035) * t849 + (t1051 * t893 + t980 * ((t883 * t981 - t975 * t990) * t964 - t863 * t981 - t963 * t1136) * t1037) * t848 + (t1049 * t892 + t978 * ((t882 * t979 - t973 * t990) * t964 - t863 * t979 - t963 * t1137) * t1036) * t847 + (t1048 * t1140 + t1050 * t1142 + t1052 * t1141) * t1012;];
taucX  = t1;
