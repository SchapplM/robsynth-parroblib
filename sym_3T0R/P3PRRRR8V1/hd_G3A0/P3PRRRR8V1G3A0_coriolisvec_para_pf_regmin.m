% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRRR8V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V1G3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:37
% EndTime: 2020-08-06 17:16:46
% DurationCPUTime: 9.45s
% Computational Cost: add. (27444->348), mult. (71070->801), div. (4896->17), fcn. (68635->22), ass. (0->376)
t967 = cos(pkin(3));
t971 = sin(qJ(3,3));
t1196 = t967 * t971;
t972 = sin(qJ(2,3));
t977 = cos(qJ(3,3));
t1183 = t972 * t977;
t978 = cos(qJ(2,3));
t930 = pkin(2) * t1183 - t978 * pkin(5);
t965 = sin(pkin(3));
t1123 = pkin(2) * t1196 + t930 * t965;
t1260 = 0.1e1 / t1123;
t983 = xDP(3);
t1187 = t967 * t983;
t1195 = t967 * t972;
t1202 = t965 * t977;
t964 = sin(pkin(6));
t1203 = t964 * t983;
t968 = legFrame(3,2);
t943 = sin(t968);
t946 = cos(t968);
t984 = xDP(2);
t985 = xDP(1);
t927 = -t943 * t984 + t946 * t985;
t966 = cos(pkin(6));
t891 = t927 * t966 - t1203;
t863 = ((-t927 * t1195 - t978 * t983) * t966 - (-t972 * t1187 + t978 * t927) * t964) * t971 - t891 * t1202;
t1220 = t863 * t1260;
t1261 = 0.1e1 / t1123 ^ 2;
t1190 = t967 * t978;
t1247 = pkin(2) * t977;
t939 = t983 * t966;
t894 = t964 * t927 + t939;
t854 = -(t891 * t1190 - t972 * t894) * t1247 - (t891 * t1195 + t978 * t894) * pkin(5);
t1231 = t854 * t1261;
t1095 = t1220 * t1231;
t1226 = t863 ^ 2 * t1261;
t1157 = t1260 * t1226;
t954 = 0.1e1 / t977 ^ 2;
t1208 = t954 * t971;
t879 = pkin(5) * (t966 * t1195 + t964 * t978) + (t966 * t1190 - t964 * t972) * t1247;
t1185 = t971 * t972;
t1029 = t967 * t1185 + t1202;
t1184 = t971 * t978;
t888 = t1029 * t966 + t964 * t1184;
t1271 = t1208 * (-0.2e1 * t888 * t1095 + t1157 * t879);
t952 = t977 ^ 2;
t953 = 0.1e1 / t977;
t955 = t953 / t952;
t1040 = (0.2e1 * t952 - 0.1e1) * t955 * t1095;
t1217 = t1260 * t953;
t1269 = t879 * t1217;
t857 = t954 * t1226;
t848 = t857 - 0.2e1 * t1226;
t1270 = 0.2e1 * t888 * t1040 + t1269 * t848;
t1216 = t1260 * t971;
t1140 = t953 * t1216;
t988 = 0.1e1 / pkin(2);
t1215 = t1260 * t988;
t1160 = t854 * t1215;
t1061 = t953 * t965 * t1160;
t1141 = t1260 * t1215;
t1094 = t863 * t1141;
t1176 = pkin(5) * t1220;
t1120 = t971 * t1176;
t1232 = t854 * t1260;
t986 = pkin(5) ^ 2;
t987 = pkin(2) ^ 2;
t839 = -pkin(5) * t854 * t1140 + (t953 * t986 + t977 * t987) * t1220;
t827 = (-t967 * t839 * t1094 - (-t971 * t930 * t1061 + t967 * (-t953 * t1120 + t977 * t1232)) * t854 * t1141) * t954;
t1077 = t827 * t1140;
t1134 = t965 * t1183;
t1135 = t965 * t1185;
t1151 = t1260 * t1220;
t1201 = t965 * t978;
t842 = (t1120 - t1232) * t953;
t824 = (-(t967 * t842 + (pkin(2) * (t967 * t1160 + t1201 * t1220) * t952 - (t854 * t1216 - t1176) * t1134) * t953) * t1151 + (-t978 * t1061 + (t1135 + (t953 - t977) * t967) * t1220) * t1231) * t954;
t1115 = t824 * t1140;
t851 = t854 ^ 2;
t1235 = t851 * t954;
t989 = 0.1e1 / pkin(2) ^ 2;
t1163 = t989 * t1235;
t899 = t1260 * t1261;
t1267 = t1115 * t879 * t988 + (t899 * t1163 + t1077) * t888;
t1101 = t851 * t899 * t955 * t971;
t1241 = t827 * t1260;
t1266 = t1215 * t824 * t879 - (t989 * t1101 - t1241) * t888;
t1244 = t824 * t1260;
t830 = (t839 * t1151 - t842 * t1231) * t953;
t1238 = t830 * t1260;
t973 = sin(qJ(3,2));
t1194 = t967 * t973;
t974 = sin(qJ(2,2));
t979 = cos(qJ(3,2));
t1180 = t974 * t979;
t980 = cos(qJ(2,2));
t931 = pkin(2) * t1180 - t980 * pkin(5);
t1122 = pkin(2) * t1194 + t931 * t965;
t1262 = 0.1e1 / t1122;
t1193 = t967 * t974;
t1200 = t965 * t979;
t969 = legFrame(2,2);
t944 = sin(t969);
t947 = cos(t969);
t928 = -t944 * t984 + t947 * t985;
t892 = t928 * t966 - t1203;
t864 = ((-t928 * t1193 - t980 * t983) * t966 - (-t974 * t1187 + t980 * t928) * t964) * t973 - t892 * t1200;
t1219 = t864 * t1262;
t1150 = t1262 * t1219;
t1263 = 0.1e1 / t1122 ^ 2;
t1189 = t967 * t980;
t1246 = pkin(2) * t979;
t895 = t964 * t928 + t939;
t855 = -(t892 * t1189 - t974 * t895) * t1246 - (t892 * t1193 + t980 * t895) * pkin(5);
t1229 = t855 * t1263;
t957 = 0.1e1 / t979;
t1207 = t957 * t973;
t1138 = t1262 * t1207;
t840 = -pkin(5) * t855 * t1138 + (t957 * t986 + t979 * t987) * t1219;
t1175 = pkin(5) * t1219;
t1119 = t973 * t1175;
t1230 = t855 * t1262;
t843 = (t1119 - t1230) * t957;
t831 = (t840 * t1150 - t843 * t1229) * t957;
t1237 = t831 * t1262;
t975 = sin(qJ(3,1));
t1192 = t967 * t975;
t976 = sin(qJ(2,1));
t981 = cos(qJ(3,1));
t1177 = t976 * t981;
t982 = cos(qJ(2,1));
t932 = pkin(2) * t1177 - t982 * pkin(5);
t1121 = pkin(2) * t1192 + t932 * t965;
t1264 = 0.1e1 / t1121;
t1191 = t967 * t976;
t1198 = t965 * t981;
t970 = legFrame(1,2);
t945 = sin(t970);
t948 = cos(t970);
t929 = -t945 * t984 + t948 * t985;
t893 = t929 * t966 - t1203;
t865 = ((-t929 * t1191 - t982 * t983) * t966 - (-t976 * t1187 + t982 * t929) * t964) * t975 - t893 * t1198;
t1218 = t865 * t1264;
t1149 = t1264 * t1218;
t1265 = 0.1e1 / t1121 ^ 2;
t1188 = t967 * t982;
t1245 = pkin(2) * t981;
t896 = t964 * t929 + t939;
t856 = -(t893 * t1188 - t976 * t896) * t1245 - (t893 * t1191 + t982 * t896) * pkin(5);
t1227 = t856 * t1265;
t961 = 0.1e1 / t981;
t1205 = t961 * t975;
t1136 = t1264 * t1205;
t841 = -pkin(5) * t856 * t1136 + (t961 * t986 + t981 * t987) * t1218;
t1174 = pkin(5) * t1218;
t1118 = t975 * t1174;
t1228 = t856 * t1264;
t844 = (t1118 - t1228) * t961;
t832 = (t841 * t1149 - t844 * t1227) * t961;
t1236 = t832 * t1264;
t1224 = t864 ^ 2 * t1263;
t1222 = t865 ^ 2 * t1265;
t1212 = t1262 * t988;
t1159 = t855 * t1212;
t1060 = t957 * t965 * t1159;
t1132 = t965 * t1180;
t1182 = t973 * t974;
t1133 = t965 * t1182;
t1199 = t965 * t980;
t1213 = t1262 * t973;
t956 = t979 ^ 2;
t958 = 0.1e1 / t979 ^ 2;
t825 = (-(t967 * t843 + (pkin(2) * (t967 * t1159 + t1199 * t1219) * t956 - (t855 * t1213 - t1175) * t1132) * t957) * t1150 - (t980 * t1060 + (-t1133 + (-t957 + t979) * t967) * t1219) * t1229) * t958;
t883 = (t964 * t1189 + t966 * t974) * t1246 + (t964 * t1193 - t966 * t980) * pkin(5);
t1259 = t825 * t883;
t1209 = t1264 * t988;
t1158 = t856 * t1209;
t1059 = t961 * t965 * t1158;
t1130 = t965 * t1177;
t1179 = t975 * t976;
t1131 = t965 * t1179;
t1197 = t965 * t982;
t1210 = t1264 * t975;
t960 = t981 ^ 2;
t962 = 0.1e1 / t981 ^ 2;
t826 = (-(t967 * t844 + (pkin(2) * (t967 * t1158 + t1197 * t1218) * t960 - (t856 * t1210 - t1174) * t1130) * t961) * t1149 + (-t982 * t1059 + (t1131 + (t961 - t981) * t967) * t1218) * t1227) * t962;
t884 = (t964 * t1188 + t966 * t976) * t1245 + (t964 * t1191 - t966 * t982) * pkin(5);
t1258 = t826 * t884;
t858 = t958 * t1224;
t859 = t962 * t1222;
t1073 = 0.2e1 * t1218 * t1227;
t1153 = t1264 * t1222;
t1204 = t962 * t975;
t881 = (t966 * t1188 - t964 * t976) * t1245 + pkin(5) * (t966 * t1191 + t964 * t982);
t1027 = t967 * t1179 + t1198;
t1178 = t975 * t982;
t890 = t1027 * t966 + t964 * t1178;
t1255 = (t890 * t1073 - t881 * t1153) * t1204;
t1074 = 0.2e1 * t1219 * t1229;
t1155 = t1262 * t1224;
t1206 = t958 * t973;
t880 = (t966 * t1189 - t964 * t974) * t1246 + pkin(5) * (t966 * t1193 + t964 * t980);
t1028 = t967 * t1182 + t1200;
t1181 = t973 * t980;
t889 = t1028 * t966 + t964 * t1181;
t1254 = (t889 * t1074 - t880 * t1155) * t1206;
t1137 = t1264 * t1209;
t1092 = t865 * t1137;
t829 = (-t967 * t841 * t1092 - (-t975 * t932 * t1059 + t967 * (-t961 * t1118 + t981 * t1228)) * t856 * t1137) * t962;
t1075 = t829 * t1136;
t1084 = t881 * t826 * t1209;
t853 = t856 ^ 2;
t1233 = t853 * t962;
t1161 = t989 * t1233;
t905 = t1264 * t1265;
t1253 = (t905 * t1161 + t1075) * t890 + t1084 * t1205;
t1139 = t1262 * t1212;
t1093 = t864 * t1139;
t828 = (-t967 * t840 * t1093 - (-t973 * t931 * t1060 + t967 * (-t957 * t1119 + t979 * t1230)) * t855 * t1139) * t958;
t1076 = t828 * t1138;
t1087 = t880 * t825 * t1212;
t852 = t855 ^ 2;
t1234 = t852 * t958;
t1162 = t989 * t1234;
t902 = t1262 * t1263;
t1252 = (t902 * t1162 + t1076) * t889 + t1087 * t1207;
t963 = t961 / t960;
t1097 = t853 * t905 * t963 * t975;
t1239 = t829 * t1264;
t1251 = (t989 * t1097 - t1239) * t890 - t1084;
t959 = t957 / t956;
t1099 = t852 * t902 * t959 * t973;
t1240 = t828 * t1262;
t1250 = (t989 * t1099 - t1240) * t889 - t1087;
t885 = t1029 * t964 - t966 * t1184;
t1249 = 0.2e1 * t885;
t1248 = pkin(2) * t965;
t1243 = t825 * t1262;
t1242 = t826 * t1264;
t1214 = t1262 * t957;
t1211 = t1264 * t961;
t1186 = t967 * t989;
t1173 = t824 * t1217;
t1172 = t824 * t1216;
t1170 = t825 * t1214;
t1169 = t825 * t1213;
t1168 = t826 * t1211;
t1167 = t826 * t1210;
t1166 = t978 * t1238;
t1165 = t980 * t1237;
t1164 = t982 * t1236;
t1156 = t1260 * t857;
t1154 = t1262 * t858;
t1152 = t1264 * t859;
t1146 = t880 * t1214;
t1145 = t881 * t1211;
t882 = (t964 * t1190 + t966 * t972) * t1247 + (t964 * t1195 - t966 * t978) * pkin(5);
t1144 = t882 * t1217;
t1143 = t883 * t1214;
t1142 = t884 * t1211;
t1026 = -0.2e1 * t854 * t978 * t1094;
t1100 = t1261 * t1163;
t821 = t824 * t1201 + t967 * t827;
t845 = t857 + t1100;
t1129 = (-t827 * t1135 + t821 * t977 + (t1026 * t1208 - t845 * t1183) * t965 - t1100 * t1196) * t1260;
t1128 = (-t827 * t1134 - t821 * t971 + (t953 * t1026 + t845 * t1185) * t965 - t953 * t851 * t1261 * t1186) * t1260;
t1025 = -0.2e1 * t855 * t980 * t1093;
t1098 = t1263 * t1162;
t822 = t825 * t1199 + t967 * t828;
t846 = t858 + t1098;
t1127 = (-t828 * t1133 + t822 * t979 + (t1025 * t1206 - t846 * t1180) * t965 - t1098 * t1194) * t1262;
t1126 = (-t828 * t1132 - t822 * t973 + (t957 * t1025 + t846 * t1182) * t965 - t957 * t852 * t1263 * t1186) * t1262;
t1023 = -0.2e1 * t856 * t982 * t1092;
t1096 = t1265 * t1161;
t823 = t826 * t1197 + t967 * t829;
t847 = t859 + t1096;
t1125 = (-t829 * t1131 + t823 * t981 + (t1023 * t1204 - t847 * t1177) * t965 - t1096 * t1192) * t1264;
t1124 = (-t829 * t1130 - t823 * t975 + (t961 * t1023 + t847 * t1179) * t965 - t961 * t853 * t1265 * t1186) * t1264;
t1117 = t888 * t1173;
t1116 = t888 * t1172;
t1114 = t889 * t1170;
t1113 = t889 * t1169;
t1112 = t890 * t1168;
t1111 = t890 * t1167;
t1110 = t830 * t1144;
t1109 = t888 * t1166;
t1108 = t830 * t972 * t1217;
t1107 = t831 * t1143;
t1106 = t889 * t1165;
t1105 = t831 * t974 * t1214;
t1104 = t832 * t1142;
t1103 = t890 * t1164;
t1102 = t832 * t976 * t1211;
t1091 = t946 * t1269;
t1090 = t943 * t1269;
t1089 = t944 * t1146;
t1088 = t947 * t1146;
t1086 = t945 * t1145;
t1085 = t948 * t1145;
t1083 = t885 * t1173;
t1082 = t885 * t1166;
t886 = t1028 * t964 - t966 * t1181;
t1081 = t886 * t1170;
t1080 = t886 * t1165;
t887 = t1027 * t964 - t966 * t1178;
t1079 = t887 * t1168;
t1078 = t887 * t1164;
t1071 = t830 * t1091;
t1070 = t830 * t1090;
t1069 = t888 * t1108;
t1068 = t831 * t1089;
t1067 = t831 * t1088;
t1066 = t889 * t1105;
t1065 = t832 * t1086;
t1064 = t832 * t1085;
t1063 = t890 * t1102;
t1058 = t953 * t1082;
t1057 = t957 * t1080;
t1056 = t961 * t1078;
t1055 = t943 * t1117;
t1054 = t943 * t1109;
t1053 = t944 * t1114;
t1052 = t944 * t1106;
t1051 = t945 * t1112;
t1050 = t945 * t1103;
t1049 = t946 * t1117;
t1048 = t946 * t1109;
t1047 = t947 * t1114;
t1046 = t947 * t1106;
t1045 = t948 * t1112;
t1044 = t948 * t1103;
t1043 = t971 * t1248 - t930 * t967;
t1042 = t973 * t1248 - t931 * t967;
t1041 = t975 * t1248 - t932 * t967;
t1037 = t953 * t1054;
t1036 = t957 * t1052;
t1035 = t961 * t1050;
t1034 = t953 * t1048;
t1033 = t957 * t1046;
t1032 = t961 * t1044;
t1024 = (0.2e1 * t956 - 0.1e1) * t959 * t1074;
t1022 = (0.2e1 * t960 - 0.1e1) * t963 * t1073;
t1020 = -t978 * t1156 - t972 * t1244;
t1019 = -t972 * t1156 + t978 * t1244;
t1018 = -t980 * t1154 - t974 * t1243;
t1017 = -t974 * t1154 + t980 * t1243;
t1016 = -t982 * t1152 - t976 * t1242;
t1015 = -t976 * t1152 + t982 * t1242;
t849 = t858 - 0.2e1 * t1224;
t1006 = t889 * t1024 + t849 * t1146;
t850 = t859 - 0.2e1 * t1222;
t1005 = t890 * t1022 + t850 * t1145;
t951 = t975 ^ 2;
t950 = t973 ^ 2;
t949 = t971 ^ 2;
t935 = pkin(5) * t976 + t982 * t1245;
t934 = pkin(5) * t974 + t980 * t1246;
t933 = pkin(5) * t972 + t978 * t1247;
t926 = t1130 + t1192;
t925 = t967 * t981 - t1131;
t924 = t1132 + t1194;
t923 = t967 * t979 - t1133;
t922 = t1134 + t1196;
t921 = t967 * t977 - t1135;
t877 = t1041 * t966 - t964 * t935;
t876 = t1042 * t966 - t964 * t934;
t875 = t1043 * t966 - t964 * t933;
t874 = t1041 * t964 + t966 * t935;
t873 = t1042 * t964 + t966 * t934;
t872 = t1043 * t964 + t966 * t933;
t871 = t1121 * t948 - t874 * t945;
t870 = t1121 * t945 + t874 * t948;
t869 = t1122 * t947 - t873 * t944;
t868 = t1122 * t944 + t873 * t947;
t867 = t1123 * t946 - t872 * t943;
t866 = t1123 * t943 + t872 * t946;
t1 = [t870 * t1236 + t868 * t1237 + t866 * t1238, -t1045 - t1047 - t1049, (t1015 * t870 + t1017 * t868 + t1019 * t866 - t1032 - t1033 - t1034) * t965, (t1016 * t870 + t1018 * t868 + t1020 * t866 + t948 * t1063 + t947 * t1066 + t946 * t1069) * t965, -t949 * t1049 - t950 * t1047 - t951 * t1045 + (-t947 * t1254 - t948 * t1255 + t946 * t1271) * t988, -0.2e1 * t946 * t1116 - 0.2e1 * t947 * t1113 - 0.2e1 * t948 * t1111 + (-t1005 * t948 - t1006 * t947 - t1270 * t946) * t988, -t1252 * t947 - t1253 * t948 - t1267 * t946, t1250 * t947 + t1251 * t948 - t1266 * t946, (-t1085 * t829 - t1088 * t828 - t1091 * t827) * t988, t870 * t1125 + t868 * t1127 + t866 * t1129 + (-t1064 * t925 - t1067 * t923 - t1071 * t921) * t988 + (-t1044 - t1046 - t1048) * t965, t870 * t1124 + t868 * t1126 + t866 * t1128 + (t1064 * t926 + t1067 * t924 + t1071 * t922) * t988 + (t1032 * t975 + t1033 * t973 + t1034 * t971) * t965, 0; t871 * t1236 + t869 * t1237 + t867 * t1238, t1051 + t1053 + t1055, (t1015 * t871 + t1017 * t869 + t1019 * t867 + t1035 + t1036 + t1037) * t965, (t1016 * t871 + t1018 * t869 + t1020 * t867 - t945 * t1063 - t944 * t1066 - t943 * t1069) * t965, t949 * t1055 + t950 * t1053 + t951 * t1051 + (t944 * t1254 + t945 * t1255 - t943 * t1271) * t988, 0.2e1 * t943 * t1116 + 0.2e1 * t944 * t1113 + 0.2e1 * t945 * t1111 + (t1005 * t945 + t1006 * t944 + t1270 * t943) * t988, t1252 * t944 + t1253 * t945 + t1267 * t943, -t1250 * t944 - t1251 * t945 + t1266 * t943, (t1086 * t829 + t1089 * t828 + t1090 * t827) * t988, t871 * t1125 + t869 * t1127 + t867 * t1129 + (t1065 * t925 + t1068 * t923 + t1070 * t921) * t988 + (t1050 + t1052 + t1054) * t965, t871 * t1124 + t869 * t1126 + t867 * t1128 + (-t1065 * t926 - t1068 * t924 - t1070 * t922) * t988 + (-t1035 * t975 - t1036 * t973 - t1037 * t971) * t965, 0; t877 * t1236 + t876 * t1237 + t875 * t1238, t1079 + t1081 + t1083, (t1015 * t877 + t1017 * t876 + t1019 * t875 + t1056 + t1057 + t1058) * t965, (t1016 * t877 + t1018 * t876 + t1020 * t875 - t887 * t1102 - t886 * t1105 - t885 * t1108) * t965, t949 * t1083 + t950 * t1081 + t951 * t1079 + ((t1073 * t887 - t1153 * t884) * t1204 + (t1074 * t886 - t1155 * t883) * t1206 + (t1095 * t1249 - t1157 * t882) * t1208) * t988, t1172 * t1249 + 0.2e1 * t886 * t1169 + 0.2e1 * t887 * t1167 + (t1022 * t887 + t1024 * t886 + t1040 * t1249 + t1142 * t850 + t1143 * t849 + t1144 * t848) * t988, t885 * t1077 + t886 * t1076 + t887 * t1075 + (t887 * t905 * t1233 + t886 * t902 * t1234 + t885 * t899 * t1235) * t989 + (t1115 * t882 + t1136 * t1258 + t1138 * t1259) * t988, t885 * t1241 + t886 * t1240 + t887 * t1239 + (-t1097 * t887 - t1099 * t886 - t1101 * t885) * t989 + (t882 * t1244 + t1258 * t1264 + t1259 * t1262) * t988, (t1142 * t829 + t1143 * t828 + t1144 * t827) * t988, t877 * t1125 + t876 * t1127 + t875 * t1129 + (t1104 * t925 + t1107 * t923 + t1110 * t921) * t988 + (t1078 + t1080 + t1082) * t965, t877 * t1124 + t876 * t1126 + t875 * t1128 + (-t1104 * t926 - t1107 * t924 - t1110 * t922) * t988 + (-t1056 * t975 - t1057 * t973 - t1058 * t971) * t965, 0;];
tau_reg  = t1;
