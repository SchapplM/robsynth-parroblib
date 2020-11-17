% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V2G1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:03:24
% EndTime: 2020-08-06 21:03:37
% DurationCPUTime: 13.10s
% Computational Cost: add. (73571->481), mult. (127527->882), div. (12054->21), fcn. (70095->70), ass. (0->407)
t1081 = pkin(5) + qJ(3,1);
t1049 = -pkin(6) - t1081;
t1114 = t1049 ^ 2;
t1078 = legFrame(1,3);
t1024 = t1078 + qJ(1,1);
t1103 = 2 * qJ(2,1);
t1062 = sin(t1103);
t1065 = cos(t1103);
t1096 = xDP(1);
t1074 = cos(pkin(7));
t1028 = t1074 * pkin(3);
t1097 = 0.2e1 * pkin(7);
t1106 = pkin(3) ^ 2;
t1342 = t1106 / 0.2e1;
t960 = cos(t1097) * t1342 + pkin(2) * (t1028 + pkin(2) / 0.2e1);
t1355 = 0.2e1 * t960;
t1228 = t1096 * t1355;
t1095 = xDP(2);
t1229 = t1095 * t1355;
t1092 = cos(qJ(2,1));
t1376 = t1028 + pkin(2);
t1275 = t1092 * t1376;
t1230 = 0.2e1 * t1275;
t1073 = sin(pkin(7));
t1086 = sin(qJ(2,1));
t1252 = t1073 * t1086;
t1256 = 0.2e1 * pkin(3);
t983 = pkin(2) * t1073 + sin(t1097) * pkin(3) / 0.2e1;
t1270 = t1096 * t983;
t1271 = t1095 * t983;
t1003 = pkin(2) * t1028;
t1108 = pkin(2) ^ 2;
t1046 = t1106 + t1108;
t1380 = 0.2e1 * t1003 + t1046;
t1306 = t1380 * t1096;
t1307 = t1380 * t1095;
t1346 = pkin(3) * t983;
t1094 = xDP(3);
t1378 = 2 * t1094;
t1336 = pkin(3) * t1073;
t943 = t1086 * t1376 + t1092 * t1336;
t1045 = pkin(1) * t1095;
t973 = -t1049 * t1096 + t1045;
t1044 = t1096 * pkin(1);
t976 = t1049 * t1095 + t1044;
t904 = (t1065 * t1228 + t976 * t1230 + (-t1062 * t1270 - t976 * t1252) * t1256 + t1306) * cos(t1024) + (t1065 * t1229 + t973 * t1230 + (-t1062 * t1271 - t973 * t1252) * t1256 + t1307) * sin(t1024) + (pkin(1) * t943 + t1062 * t960 + t1065 * t1346) * t1378;
t1225 = pkin(3) * t1252;
t940 = 0.1e1 / (-t1225 + t1275);
t1385 = 0.1e1 / t1114 * t904 * t940;
t1080 = pkin(5) + qJ(3,2);
t1048 = -pkin(6) - t1080;
t1113 = t1048 ^ 2;
t1077 = legFrame(2,3);
t1023 = t1077 + qJ(1,2);
t1101 = 2 * qJ(2,2);
t1061 = sin(t1101);
t1064 = cos(t1101);
t1090 = cos(qJ(2,2));
t1278 = t1090 * t1376;
t1231 = 0.2e1 * t1278;
t1084 = sin(qJ(2,2));
t1253 = t1073 * t1084;
t942 = t1084 * t1376 + t1090 * t1336;
t972 = -t1048 * t1096 + t1045;
t975 = t1048 * t1095 + t1044;
t903 = (t1064 * t1228 + t975 * t1231 + (-t1061 * t1270 - t975 * t1253) * t1256 + t1306) * cos(t1023) + (t1064 * t1229 + t972 * t1231 + (-t1061 * t1271 - t972 * t1253) * t1256 + t1307) * sin(t1023) + (pkin(1) * t942 + t1061 * t960 + t1064 * t1346) * t1378;
t1226 = pkin(3) * t1253;
t939 = 0.1e1 / (-t1226 + t1278);
t1384 = 0.1e1 / t1113 * t903 * t939;
t1079 = pkin(5) + qJ(3,3);
t1047 = -pkin(6) - t1079;
t1112 = t1047 ^ 2;
t1076 = legFrame(3,3);
t1022 = t1076 + qJ(1,3);
t1099 = 2 * qJ(2,3);
t1060 = sin(t1099);
t1063 = cos(t1099);
t1088 = cos(qJ(2,3));
t1281 = t1088 * t1376;
t1232 = 0.2e1 * t1281;
t1082 = sin(qJ(2,3));
t1254 = t1073 * t1082;
t941 = t1082 * t1376 + t1088 * t1336;
t971 = -t1047 * t1096 + t1045;
t974 = t1047 * t1095 + t1044;
t902 = (t1063 * t1228 + t974 * t1232 + (-t1060 * t1270 - t974 * t1254) * t1256 + t1306) * cos(t1022) + (t1063 * t1229 + t971 * t1232 + (-t1060 * t1271 - t971 * t1254) * t1256 + t1307) * sin(t1022) + (pkin(1) * t941 + t1060 * t960 + t1063 * t1346) * t1378;
t1227 = pkin(3) * t1254;
t938 = 0.1e1 / (-t1227 + t1281);
t1383 = 0.1e1 / t1112 * t902 * t938;
t1382 = 0.3e1 * t1108 + 0.6e1 * t1106;
t1381 = 0.6e1 * t1108 + 0.3e1 * t1106;
t1041 = t1088 * pkin(2);
t1055 = qJ(2,3) + pkin(7);
t1017 = cos(t1055);
t1339 = pkin(3) * t1017;
t980 = t1041 + t1339;
t1042 = t1090 * pkin(2);
t1057 = qJ(2,2) + pkin(7);
t1019 = cos(t1057);
t1338 = pkin(3) * t1019;
t981 = t1042 + t1338;
t1043 = t1092 * pkin(2);
t1059 = qJ(2,1) + pkin(7);
t1021 = cos(t1059);
t1337 = pkin(3) * t1021;
t982 = t1043 + t1337;
t1379 = 0.2e1 * pkin(1);
t1377 = -pkin(5) / 0.2e1;
t999 = t1041 + pkin(1);
t1000 = t1042 + pkin(1);
t1001 = t1043 + pkin(1);
t1035 = 0.1e1 / t1047;
t1324 = t938 * t941;
t1029 = sin(t1076);
t1032 = cos(t1076);
t1083 = sin(qJ(1,3));
t1089 = cos(qJ(1,3));
t953 = -t1029 * t1083 + t1032 * t1089;
t956 = t1029 * t1089 + t1032 * t1083;
t923 = (t1094 * t1324 + t1095 * t956 + t1096 * t953) * t1035;
t962 = 0.1e1 / t980;
t1327 = t923 * t962;
t1037 = 0.1e1 / t1048;
t1323 = t939 * t942;
t1030 = sin(t1077);
t1033 = cos(t1077);
t1085 = sin(qJ(1,2));
t1091 = cos(qJ(1,2));
t954 = -t1030 * t1085 + t1033 * t1091;
t957 = t1030 * t1091 + t1033 * t1085;
t924 = (t1094 * t1323 + t1095 * t957 + t1096 * t954) * t1037;
t965 = 0.1e1 / t981;
t1326 = t924 * t965;
t1039 = 0.1e1 / t1049;
t1322 = t940 * t943;
t1031 = sin(t1078);
t1034 = cos(t1078);
t1087 = sin(qJ(1,1));
t1093 = cos(qJ(1,1));
t955 = -t1031 * t1087 + t1034 * t1093;
t958 = t1031 * t1093 + t1034 * t1087;
t925 = (t1094 * t1322 + t1095 * t958 + t1096 * t955) * t1039;
t968 = 0.1e1 / t982;
t1325 = t925 * t968;
t1069 = t1088 ^ 2;
t1054 = t1099 + pkin(7);
t1016 = cos(t1054);
t1243 = t1088 * t1379;
t1318 = pkin(3) * t1379;
t996 = 0.2e1 * t1055;
t1368 = t1108 * t1063 + t1106 * cos(t996);
t1143 = (t1017 * t1318 + t1046 + (t1243 + (t1016 + t1074) * t1256) * pkin(2) + t1368) * t1327 * t1383;
t1072 = t1094 ^ 2;
t1115 = t980 ^ 2;
t963 = 0.1e1 / t1115;
t977 = t1082 * pkin(2) + pkin(3) * sin(t1055);
t1321 = t962 * t963 * t977;
t1209 = t1082 * t1321;
t1152 = t1072 * t1209;
t1146 = pkin(2) * t1152;
t1260 = t962 * t1082;
t1188 = t1094 * t1260;
t1164 = pkin(2) * t1188;
t1274 = t1094 * t962;
t1185 = t1274 / 0.2e1;
t1010 = cos(t1097 + qJ(2,3));
t1013 = cos(-pkin(7) + qJ(2,3));
t1350 = pkin(2) * pkin(3);
t1236 = sin(t1054) * t1350;
t1365 = t1108 * t1060 + t1106 * sin(t996);
t1128 = 0.2e1 * t1236 + t1365;
t1105 = pkin(3) * t1106;
t1335 = pkin(3) * t1108;
t1171 = -0.2e1 * t1105 - 0.4e1 * t1335;
t1237 = -0.2e1 * t1335;
t1340 = pkin(2) * t1106;
t1239 = -0.2e1 * t1340;
t1107 = pkin(2) * t1108;
t1353 = -0.2e1 * t1107 - 0.4e1 * t1340;
t1245 = t1003 + t1108 / 0.2e1;
t1354 = -0.4e1 * pkin(1) * (t1342 + t1245);
t1200 = (t1128 * t923 * t1047 + (t1010 * t1239 + t1013 * t1237 + t1171 * t1017 + t1088 * t1353 + t1354) * t1274) * t963 * t1094;
t1302 = t1035 * t962;
t1221 = t923 * t1302;
t1098 = 3 * qJ(2,3);
t1109 = pkin(1) ^ 2;
t1305 = t1035 * t938;
t1197 = t902 * t1305;
t1238 = -0.3e1 * t1335;
t1240 = -0.3e1 * t1340;
t1317 = -0.8e1 * t1350;
t1372 = -0.8e1 * t1003 - 0.4e1 * t1046;
t899 = t1197 / 0.2e1;
t918 = pkin(1) * t923;
t857 = (-0.4e1 * t1236 - 0.2e1 * t1365) * t1047 * t1274 - 0.4e1 * t980 * (pkin(1) * t899 - (t1109 + t1112) * t923) + (t1105 * cos(0.3e1 * t1055) + t1339 * t1381 + t1107 * cos(t1098) + t1041 * t1382 - (cos(t1097 + t1098) + t1010) * t1240 - (cos(t1098 + pkin(7)) + t1013) * t1238) * t923 + (t1016 * t1317 - 0.4e1 * t1368 + t1372) * (-t918 + t1197 / 0.4e1);
t1224 = -t857 * t1221 / 0.4e1 + t1035 * t1200 / 0.2e1 + t1143 / 0.4e1;
t1244 = pkin(5) ^ 2 + t1109;
t1286 = t1082 * t923;
t1290 = t1072 * t1380;
t1349 = t923 / 0.2e1;
t1053 = t1074 ^ 2;
t1356 = -0.2e1 * (t1053 - 0.1e1 / 0.2e1) * t1106 - 0.2e1 * t1245;
t1371 = (t1053 - 0.1e1) * t1106;
t893 = -t918 + t899;
t869 = -t963 * t1035 / (t1041 + (t1074 * t1088 - t1254) * pkin(3)) * t1290 + t1349 * t1383 - (t893 * t1227 - (t1069 * t1356 + t1371) * t923 + (-0.2e1 * t923 * t1227 - t893) * t1281) * t1305 * t923;
t1345 = t869 * pkin(1);
t1363 = 0.2e1 * pkin(5);
t1316 = 0.2e1 * (t1345 - t1143 / 0.8e1 - (-t857 * t1327 / 0.8e1 + t1200 / 0.4e1) * t1035) * t1041 - t1079 * t1146 - pkin(1) * t1224 + (t1069 * t1108 + (t1363 + qJ(3,3)) * qJ(3,3) + t1244) * t869 - 0.2e1 * (-pkin(2) * t1286 + t1079 * t1185) * t1274 * t1041 + 0.2e1 * t923 * (pkin(1) * t1164 + t1079 * t899);
t1070 = t1090 ^ 2;
t1056 = t1101 + pkin(7);
t1018 = cos(t1056);
t1242 = t1090 * t1379;
t997 = 0.2e1 * t1057;
t1369 = t1108 * t1064 + t1106 * cos(t997);
t1142 = (t1019 * t1318 + t1046 + (t1242 + (t1018 + t1074) * t1256) * pkin(2) + t1369) * t1326 * t1384;
t1118 = t981 ^ 2;
t966 = 0.1e1 / t1118;
t978 = t1084 * pkin(2) + pkin(3) * sin(t1057);
t1320 = t965 * t966 * t978;
t1208 = t1084 * t1320;
t1151 = t1072 * t1208;
t1145 = pkin(2) * t1151;
t1259 = t965 * t1084;
t1187 = t1094 * t1259;
t1163 = pkin(2) * t1187;
t1273 = t1094 * t965;
t1184 = t1273 / 0.2e1;
t1011 = cos(t1097 + qJ(2,2));
t1014 = cos(-pkin(7) + qJ(2,2));
t1235 = sin(t1056) * t1350;
t1366 = t1108 * t1061 + t1106 * sin(t997);
t1127 = 0.2e1 * t1235 + t1366;
t1199 = (t1127 * t924 * t1048 + (t1011 * t1239 + t1014 * t1237 + t1171 * t1019 + t1090 * t1353 + t1354) * t1273) * t966 * t1094;
t1297 = t1037 * t965;
t1217 = t924 * t1297;
t1100 = 3 * qJ(2,2);
t1300 = t1037 * t939;
t1196 = t903 * t1300;
t900 = t1196 / 0.2e1;
t919 = pkin(1) * t924;
t858 = (-0.4e1 * t1235 - 0.2e1 * t1366) * t1048 * t1273 - 0.4e1 * t981 * (pkin(1) * t900 - (t1109 + t1113) * t924) + (t1105 * cos(0.3e1 * t1057) + t1338 * t1381 + t1107 * cos(t1100) + t1042 * t1382 - (cos(t1097 + t1100) + t1011) * t1240 - (cos(t1100 + pkin(7)) + t1014) * t1238) * t924 + (t1018 * t1317 - 0.4e1 * t1369 + t1372) * (-t919 + t1196 / 0.4e1);
t1223 = -t858 * t1217 / 0.4e1 + t1037 * t1199 / 0.2e1 + t1142 / 0.4e1;
t1285 = t1084 * t924;
t1348 = t924 / 0.2e1;
t895 = -t919 + t900;
t870 = -t966 * t1037 / (t1042 + (t1074 * t1090 - t1253) * pkin(3)) * t1290 + t1348 * t1384 - (t895 * t1226 - (t1070 * t1356 + t1371) * t924 + (-0.2e1 * t924 * t1226 - t895) * t1278) * t1300 * t924;
t1344 = t870 * pkin(1);
t1315 = 0.2e1 * (t1344 - t1142 / 0.8e1 - (-t858 * t1326 / 0.8e1 + t1199 / 0.4e1) * t1037) * t1042 - t1080 * t1145 - pkin(1) * t1223 + (t1070 * t1108 + (t1363 + qJ(3,2)) * qJ(3,2) + t1244) * t870 - 0.2e1 * (-pkin(2) * t1285 + t1080 * t1184) * t1273 * t1042 + 0.2e1 * t924 * (pkin(1) * t1163 + t1080 * t900);
t1071 = t1092 ^ 2;
t1058 = t1103 + pkin(7);
t1020 = cos(t1058);
t1241 = t1092 * t1379;
t998 = 0.2e1 * t1059;
t1370 = t1108 * t1065 + t1106 * cos(t998);
t1141 = (t1021 * t1318 + t1046 + (t1241 + (t1020 + t1074) * t1256) * pkin(2) + t1370) * t1325 * t1385;
t1121 = t982 ^ 2;
t969 = 0.1e1 / t1121;
t979 = t1086 * pkin(2) + pkin(3) * sin(t1059);
t1319 = t968 * t969 * t979;
t1207 = t1086 * t1319;
t1150 = t1072 * t1207;
t1144 = pkin(2) * t1150;
t1258 = t968 * t1086;
t1186 = t1094 * t1258;
t1162 = pkin(2) * t1186;
t1272 = t1094 * t968;
t1183 = t1272 / 0.2e1;
t1012 = cos(t1097 + qJ(2,1));
t1015 = cos(-pkin(7) + qJ(2,1));
t1234 = sin(t1058) * t1350;
t1367 = t1108 * t1062 + t1106 * sin(t998);
t1126 = 0.2e1 * t1234 + t1367;
t1198 = (t1126 * t925 * t1049 + (t1012 * t1239 + t1015 * t1237 + t1171 * t1021 + t1092 * t1353 + t1354) * t1272) * t969 * t1094;
t1292 = t1039 * t968;
t1213 = t925 * t1292;
t1102 = 3 * qJ(2,1);
t1295 = t1039 * t940;
t1195 = t904 * t1295;
t901 = t1195 / 0.2e1;
t917 = t925 * pkin(1);
t859 = (-0.4e1 * t1234 - 0.2e1 * t1367) * t1049 * t1272 - 0.4e1 * t982 * (pkin(1) * t901 - (t1109 + t1114) * t925) + (t1105 * cos(0.3e1 * t1059) + t1337 * t1381 + t1107 * cos(t1102) + t1043 * t1382 - (cos(t1102 + t1097) + t1012) * t1240 - (cos(t1102 + pkin(7)) + t1015) * t1238) * t925 + (t1020 * t1317 - 0.4e1 * t1370 + t1372) * (-t917 + t1195 / 0.4e1);
t1222 = -t859 * t1213 / 0.4e1 + t1039 * t1198 / 0.2e1 + t1141 / 0.4e1;
t1284 = t1086 * t925;
t1347 = t925 / 0.2e1;
t891 = -t917 + t901;
t871 = -t969 * t1039 / (t1043 + (t1074 * t1092 - t1252) * pkin(3)) * t1290 + t1347 * t1385 - (t891 * t1225 - (t1071 * t1356 + t1371) * t925 + (-0.2e1 * t925 * t1225 - t891) * t1275) * t1295 * t925;
t1343 = t871 * pkin(1);
t1314 = 0.2e1 * (t1343 - t1141 / 0.8e1 - (-t859 * t1325 / 0.8e1 + t1198 / 0.4e1) * t1039) * t1043 - t1081 * t1144 - pkin(1) * t1222 + (t1071 * t1108 + (t1363 + qJ(3,1)) * qJ(3,1) + t1244) * t871 - 0.2e1 * (-pkin(2) * t1284 + t1081 * t1183) * t1272 * t1043 + 0.2e1 * t925 * (pkin(1) * t1162 + t1081 * t901);
t1025 = 0.2e1 * t1069 - 0.1e1;
t1026 = 0.2e1 * t1070 - 0.1e1;
t1027 = 0.2e1 * t1071 - 0.1e1;
t1341 = pkin(2) * t1072;
t1334 = pkin(5) * t1072;
t1333 = (t1079 * t1349 + t1164) * t923;
t1332 = (t1080 * t1348 + t1163) * t924;
t1331 = (t1081 * t1347 + t1162) * t925;
t920 = t923 ^ 2;
t1330 = t920 * t962;
t921 = t924 ^ 2;
t1329 = t921 * t965;
t922 = t925 ^ 2;
t1328 = t922 * t968;
t1282 = t1088 * t963;
t1289 = t1079 * t869;
t1313 = t923 * t1197 - t1282 * t1341 - t1146 + 0.2e1 * t1289;
t1279 = t1090 * t966;
t1288 = t1080 * t870;
t1312 = t924 * t1196 - t1279 * t1341 - t1145 + 0.2e1 * t1288;
t1276 = t1092 * t969;
t1287 = t1081 * t871;
t1311 = t925 * t1195 - t1276 * t1341 - t1144 + 0.2e1 * t1287;
t1304 = t1035 * t953;
t1303 = t1035 * t956;
t1299 = t1037 * t954;
t1298 = t1037 * t957;
t1294 = t1039 * t955;
t1293 = t1039 * t958;
t1283 = t1088 * t962;
t1280 = t1090 * t965;
t1277 = t1092 * t968;
t1263 = t923 * t1088;
t1262 = t924 * t1090;
t1261 = t925 * t1092;
t1220 = t941 * t1305;
t1219 = t953 * t1302;
t1218 = t956 * t1302;
t1216 = t942 * t1300;
t1215 = t954 * t1297;
t1214 = t957 * t1297;
t1212 = t943 * t1295;
t1211 = t955 * t1292;
t1210 = t958 * t1292;
t1206 = t920 * t1283;
t1205 = t1088 * t1321;
t1204 = t921 * t1280;
t1203 = t1090 * t1320;
t1202 = t922 * t1277;
t1201 = t1092 * t1319;
t1066 = t1082 ^ 2;
t1191 = t1035 * t1066 * t869;
t1067 = t1084 ^ 2;
t1190 = t1037 * t1067 * t870;
t1068 = t1086 ^ 2;
t1189 = t1039 * t1068 * t871;
t1182 = t1035 * t1082 * t1088;
t1181 = t1037 * t1084 * t1090;
t1180 = t1039 * t1086 * t1092;
t851 = -t999 * t869 + t1224;
t1170 = t851 - 0.2e1 * t1333;
t852 = -t1000 * t870 + t1223;
t1169 = t852 - 0.2e1 * t1332;
t853 = -t1001 * t871 + t1222;
t1168 = t853 - 0.2e1 * t1331;
t1161 = t869 * t1220;
t1160 = t962 * t1220;
t1159 = t870 * t1216;
t1158 = t965 * t1216;
t1157 = t871 * t1212;
t1156 = t968 * t1212;
t1155 = t1025 * t1221;
t1154 = t1026 * t1217;
t1153 = t1027 * t1213;
t1149 = t869 * t1182;
t1148 = t870 * t1181;
t1147 = t871 * t1180;
t1140 = t1182 * t1327;
t1139 = t1182 * t1324;
t1138 = t1181 * t1326;
t1137 = t1181 * t1323;
t1136 = t1180 * t1325;
t1135 = t1180 * t1322;
t1134 = t1035 * (-t1082 * t963 + t1205);
t1133 = t1035 * (t1209 + t1282);
t1132 = t1037 * (-t1084 * t966 + t1203);
t1131 = t1037 * (t1208 + t1279);
t1130 = t1039 * (-t1086 * t969 + t1201);
t1129 = t1039 * (t1207 + t1276);
t1125 = t871 * t1258 + t870 * t1259 + t869 * t1260;
t1124 = t871 * t1277 + t870 * t1280 + t869 * t1283;
t952 = t1000 * t1091 - t1048 * t1085;
t951 = -t1047 * t1083 + t1089 * t999;
t950 = t1001 * t1093 - t1049 * t1087;
t949 = t1001 * t1087 + t1049 * t1093;
t948 = t1000 * t1085 + t1048 * t1091;
t947 = t1047 * t1089 + t1083 * t999;
t937 = t1379 * t979 + t1126;
t936 = t1379 * t978 + t1127;
t935 = t1379 * t977 + t1128;
t931 = t1031 * t950 + t1034 * t949 + t958 * t1337;
t930 = t1030 * t952 + t1033 * t948 + t957 * t1338;
t929 = t1029 * t951 + t1032 * t947 + t956 * t1339;
t928 = -t1031 * t949 + t1034 * t950 + t955 * t1337;
t927 = -t1030 * t948 + t1033 * t952 + t954 * t1338;
t926 = -t1029 * t947 + t1032 * t951 + t953 * t1339;
t910 = -pkin(1) * t1261 + t1186 * t1377;
t909 = -pkin(1) * t1262 + t1187 * t1377;
t908 = -pkin(1) * t1263 + t1188 * t1377;
t907 = t1092 * pkin(5) * t1183 - pkin(1) * t1284;
t906 = t1090 * pkin(5) * t1184 - pkin(1) * t1285;
t905 = t1088 * pkin(5) * t1185 - pkin(1) * t1286;
t865 = -0.2e1 * t1086 * t1343 - t1201 * t1334;
t864 = -0.2e1 * t1084 * t1344 - t1203 * t1334;
t863 = -0.2e1 * t1082 * t1345 - t1205 * t1334;
t862 = -pkin(5) * t1150 + t871 * t1241;
t861 = -pkin(5) * t1151 + t870 * t1242;
t860 = -pkin(5) * t1152 + t869 * t1243;
t1 = [-t871 * t1294 - t870 * t1299 - t869 * t1304, 0, 0, -t953 * t1191 - t954 * t1190 - t955 * t1189 + (t1136 * t955 + t1138 * t954 + t1140 * t953) * t1378, -0.2e1 * t953 * t1149 - 0.2e1 * t954 * t1148 - 0.2e1 * t955 * t1147 + 0.2e1 * (t1153 * t955 + t1154 * t954 + t1155 * t953) * t1094, (-t1129 * t955 - t1131 * t954 - t1133 * t953) * t1072, (-t1130 * t955 - t1132 * t954 - t1134 * t953) * t1072, 0, -t860 * t1304 - t861 * t1299 - t862 * t1294 + (t907 * t1211 + t906 * t1215 + t905 * t1219) * t1378, -t863 * t1304 - t864 * t1299 - t865 * t1294 + (t910 * t1211 + t909 * t1215 + t908 * t1219) * t1378, -(t1311 * t955 - t922 * t928) * t1039 - (t1312 * t954 - t921 * t927) * t1037 - (t1313 * t953 - t920 * t926) * t1035, -(t1168 * t928 + t1314 * t955) * t1039 - (t1169 * t927 + t1315 * t954) * t1037 - (t1170 * t926 + t1316 * t953) * t1035, 0; -t871 * t1293 - t870 * t1298 - t869 * t1303, 0, 0, -t956 * t1191 - t957 * t1190 - t958 * t1189 + (t1136 * t958 + t1138 * t957 + t1140 * t956) * t1378, -0.2e1 * t956 * t1149 - 0.2e1 * t957 * t1148 - 0.2e1 * t958 * t1147 + 0.2e1 * (t1153 * t958 + t1154 * t957 + t1155 * t956) * t1094, (-t1129 * t958 - t1131 * t957 - t1133 * t956) * t1072, (-t1130 * t958 - t1132 * t957 - t1134 * t956) * t1072, 0, -t860 * t1303 - t861 * t1298 - t862 * t1293 + (t907 * t1210 + t906 * t1214 + t905 * t1218) * t1378, -t863 * t1303 - t864 * t1298 - t865 * t1293 + (t1210 * t910 + t1214 * t909 + t1218 * t908) * t1378, -(t1311 * t958 - t922 * t931) * t1039 - (t1312 * t957 - t921 * t930) * t1037 - (t1313 * t956 - t920 * t929) * t1035, -(t1168 * t931 + t1314 * t958) * t1039 - (t1169 * t930 + t1315 * t957) * t1037 - (t1170 * t929 + t1316 * t956) * t1035, 0; -t1157 - t1159 - t1161, 0, 0, -t1066 * t1161 - t1067 * t1159 - t1068 * t1157 - t1082 * t1206 - t1084 * t1204 - t1086 * t1202 + (t1135 * t1325 + t1137 * t1326 + t1139 * t1327) * t1378, -t1027 * t1328 - t1026 * t1329 - t1025 * t1330 - 0.2e1 * t871 * t1135 - 0.2e1 * t870 * t1137 - 0.2e1 * t869 * t1139 + (t1153 * t1322 + t1154 * t1323 + t1155 * t1324) * t1378, (-t1129 * t1322 - t1131 * t1323 - t1133 * t1324) * t1072 + t1125, (-t1130 * t1322 - t1132 * t1323 - t1134 * t1324) * t1072 + t1124, (0.1e1 / t1121 ^ 2 * t979 + 0.1e1 / t1118 ^ 2 * t978 + 0.1e1 / t1115 ^ 2 * t977) * t1072, -t860 * t1220 - t861 * t1216 - t862 * t1212 + (t1156 * t907 + t1158 * t906 + t1160 * t905) * t1378 - t1125 * pkin(5) + (t922 * t1258 + t921 * t1259 + t920 * t1260) * pkin(1), -t863 * t1220 - t864 * t1216 - t865 * t1212 + (t1156 * t910 + t1158 * t909 + t1160 * t908) * t1378 - t1124 * pkin(5) + (t1202 + t1204 + t1206) * pkin(1), -(-t937 * t1328 / 0.2e1 + t1311 * t1322) * t1039 - (-t936 * t1329 / 0.2e1 + t1312 * t1323) * t1037 - (-t935 * t1330 / 0.2e1 + t1313 * t1324) * t1035 - t1125 * pkin(2), ((t1319 * t1341 + (-t1287 - (-pkin(2) * t1261 + t1195 - t917) * t925) * t1086) * pkin(2) - (-t1331 + t853 / 0.2e1) * t1039 * t937) * t968 + ((t1320 * t1341 + (-t1288 - (-pkin(2) * t1262 + t1196 - t919) * t924) * t1084) * pkin(2) - (-t1332 + t852 / 0.2e1) * t1037 * t936) * t965 + ((t1321 * t1341 + (-t1289 - (-pkin(2) * t1263 + t1197 - t918) * t923) * t1082) * pkin(2) - (-t1333 + t851 / 0.2e1) * t1035 * t935) * t962 - t1316 * t1220 - t1315 * t1216 - t1314 * t1212, 0;];
tau_reg  = t1;