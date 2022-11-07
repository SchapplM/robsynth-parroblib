% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRPRR8V2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
% tau_reg [3*3x15]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V2G2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:06:40
% EndTime: 2022-11-07 13:06:51
% DurationCPUTime: 11.80s
% Computational Cost: add. (14251->609), mult. (26403->1270), div. (2274->18), fcn. (20424->23), ass. (0->552)
t1390 = 2 * pkin(1);
t1029 = cos(pkin(7));
t1389 = 0.2e1 * t1029;
t1030 = (qJ(3,3) + pkin(5));
t1388 = 2 * t1030;
t1031 = (qJ(3,2) + pkin(5));
t1387 = 2 * t1031;
t1032 = (qJ(3,1) + pkin(5));
t1386 = 2 * t1032;
t1036 = sin(qJ(2,3));
t1037 = sin(qJ(1,3));
t1042 = cos(qJ(2,3));
t1033 = legFrame(3,2);
t1002 = sin(t1033);
t1028 = sin(pkin(7));
t1357 = t1028 * pkin(3);
t1271 = t1002 * t1357;
t1005 = cos(t1033);
t1274 = t1005 * t1357;
t1305 = t1002 * t1037;
t1356 = t1029 * pkin(3);
t991 = pkin(2) + t1356;
t1312 = t991 * t1005;
t930 = (-t991 * t1305 + t1274) * t1042 + (t1037 * t1271 + t1312) * t1036;
t1302 = t1005 * t1037;
t1336 = t1002 * t991;
t933 = (t991 * t1302 + t1271) * t1042 + t1036 * (-t1037 * t1274 + t1336);
t1385 = t930 * t933;
t1038 = sin(qJ(2,2));
t1039 = sin(qJ(1,2));
t1044 = cos(qJ(2,2));
t1034 = legFrame(2,2);
t1003 = sin(t1034);
t1270 = t1003 * t1357;
t1006 = cos(t1034);
t1273 = t1006 * t1357;
t1304 = t1003 * t1039;
t1311 = t991 * t1006;
t931 = (-t991 * t1304 + t1273) * t1044 + (t1039 * t1270 + t1311) * t1038;
t1301 = t1006 * t1039;
t1334 = t1003 * t991;
t934 = (t991 * t1301 + t1270) * t1044 + t1038 * (-t1039 * t1273 + t1334);
t1384 = t931 * t934;
t1040 = sin(qJ(2,1));
t1041 = sin(qJ(1,1));
t1046 = cos(qJ(2,1));
t1035 = legFrame(1,2);
t1004 = sin(t1035);
t1269 = t1004 * t1357;
t1007 = cos(t1035);
t1272 = t1007 * t1357;
t1303 = t1004 * t1041;
t1310 = t991 * t1007;
t932 = (-t991 * t1303 + t1272) * t1046 + (t1041 * t1269 + t1310) * t1040;
t1300 = t1007 * t1041;
t1332 = t1004 * t991;
t935 = (t991 * t1300 + t1269) * t1046 + t1040 * (-t1041 * t1272 + t1332);
t1383 = t932 * t935;
t1381 = 2 * pkin(5);
t1022 = t1042 ^ 2;
t1380 = 0.2e1 * t1022;
t1024 = t1044 ^ 2;
t1379 = 0.2e1 * t1024;
t1026 = t1046 ^ 2;
t1378 = 0.2e1 * t1026;
t1377 = -0.2e1 * t1036;
t1376 = -0.2e1 * t1038;
t1375 = -0.2e1 * t1040;
t1374 = 0.2e1 * t1042;
t1373 = 0.2e1 * t1044;
t1372 = 0.2e1 * t1046;
t1371 = pkin(1) * t930;
t1370 = pkin(1) * t931;
t1369 = pkin(1) * t932;
t1368 = pkin(1) * t933;
t1367 = pkin(1) * t934;
t1366 = pkin(1) * t935;
t1015 = pkin(6) + t1030;
t1043 = cos(qJ(1,3));
t1164 = pkin(1) * t1037 - t1043 * t1015;
t1018 = t1029 ^ 2;
t1276 = pkin(3) * (t1018 - 0.1e1);
t1298 = t1028 * t1036;
t1365 = pkin(3) * (t1037 * t1276 + t1164 * t1298);
t1016 = pkin(6) + t1031;
t1045 = cos(qJ(1,2));
t1163 = pkin(1) * t1039 - t1045 * t1016;
t1297 = t1028 * t1038;
t1364 = pkin(3) * (t1039 * t1276 + t1163 * t1297);
t1017 = pkin(6) + t1032;
t1047 = cos(qJ(1,1));
t1162 = pkin(1) * t1041 - t1047 * t1017;
t1296 = t1028 * t1040;
t1363 = pkin(3) * (t1041 * t1276 + t1162 * t1296);
t1362 = t1029 / 0.2e1;
t1361 = pkin(1) * t1043;
t1360 = pkin(1) * t1045;
t1359 = pkin(1) * t1047;
t1358 = pkin(2) * t1029;
t1008 = t1042 * pkin(2);
t1009 = t1044 * pkin(2);
t1010 = t1046 * pkin(2);
t996 = 1 / t1015;
t1355 = t930 * t996;
t998 = 1 / t1016;
t1354 = t931 * t998;
t1353 = t933 * t996;
t1352 = t934 * t998;
t1268 = pkin(3) * t1298;
t966 = t991 * t1042 - t1268;
t960 = 0.1e1 / t966;
t1351 = t960 * t996;
t961 = 0.1e1 / t966 ^ 2;
t997 = 1 / t1015 ^ 2;
t1350 = t961 * t997;
t1267 = pkin(3) * t1297;
t967 = t991 * t1044 - t1267;
t962 = 0.1e1 / t967;
t1349 = t962 * t998;
t963 = 0.1e1 / t967 ^ 2;
t999 = 1 / t1016 ^ 2;
t1348 = t963 * t999;
t970 = -t1029 * t1042 + t1298;
t1347 = t970 * t997;
t1069 = -t1029 * t1044 + t1297;
t1346 = t1069 * t999;
t1295 = t1029 * t1036;
t973 = t1028 * t1042 + t1295;
t1345 = t973 * t997;
t1294 = t1029 * t1038;
t974 = t1028 * t1044 + t1294;
t1344 = t974 * t999;
t1000 = 1 / t1017;
t1343 = t1000 * t932;
t1342 = t1000 * t935;
t1266 = pkin(3) * t1296;
t968 = t991 * t1046 - t1266;
t964 = 0.1e1 / t968;
t1341 = t1000 * t964;
t1001 = 1 / t1017 ^ 2;
t965 = 0.1e1 / t968 ^ 2;
t1340 = t1001 * t965;
t1068 = -t1029 * t1046 + t1296;
t1339 = t1001 * t1068;
t1293 = t1029 * t1040;
t975 = t1028 * t1046 + t1293;
t1338 = t1001 * t975;
t982 = pkin(3) * cos(qJ(2,3) + pkin(7)) + t1008;
t976 = 0.1e1 / t982;
t1337 = t1002 * t976;
t983 = pkin(3) * cos(qJ(2,2) + pkin(7)) + t1009;
t978 = 0.1e1 / t983;
t1335 = t1003 * t978;
t984 = pkin(3) * cos(qJ(2,1) + pkin(7)) + t1010;
t980 = 0.1e1 / t984;
t1333 = t1004 * t980;
t1331 = t1005 * t976;
t1330 = t1006 * t978;
t1329 = t1007 * t980;
t1023 = t1043 ^ 2;
t1328 = t1023 * t997;
t1025 = t1045 ^ 2;
t1327 = t1025 * t999;
t1326 = t1030 * t976;
t1325 = t1031 * t978;
t1324 = t1032 * t980;
t1323 = t1036 * t976;
t1322 = t1038 * t978;
t1321 = t1040 * t980;
t1320 = t1042 * t996;
t1319 = t1042 * t997;
t1318 = t1043 * t996;
t1317 = t1043 * t997;
t1316 = t1044 * t998;
t1315 = t1044 * t999;
t1314 = t1045 * t998;
t1313 = t1045 * t999;
t1308 = t1000 * t1046;
t1307 = t1000 * t1047;
t1306 = t1001 * t1047;
t1027 = t1047 ^ 2;
t1299 = t1027 * t1001;
t1292 = pkin(1) ^ 2 + pkin(5) ^ 2;
t1291 = pkin(2) * t1356;
t1290 = pkin(2) * t1022 * t996;
t1289 = pkin(2) * t1024 * t998;
t1288 = pkin(5) * t1042 * t976;
t1287 = pkin(5) * t1044 * t978;
t1286 = pkin(5) * t1046 * t980;
t1285 = pkin(2) * t1337;
t1284 = pkin(2) * t1335;
t1283 = pkin(2) * t1333;
t1282 = pkin(2) * t1331;
t1281 = pkin(2) * t1330;
t1280 = pkin(2) * t1329;
t1279 = -0.2e1 * t1298;
t1278 = -0.2e1 * t1297;
t1277 = -0.2e1 * t1296;
t1275 = pkin(2) * t1000 * t1026;
t1265 = t930 * t1351;
t1264 = t931 * t1349;
t1263 = t933 * t1351;
t1262 = t934 * t1349;
t954 = t1037 * t1015 + (pkin(1) + t982) * t1043;
t1261 = t954 * t1347;
t1260 = t954 * t1345;
t955 = t1039 * t1016 + (pkin(1) + t983) * t1045;
t1259 = t955 * t1346;
t1258 = t955 * t1344;
t1257 = t961 * t1347;
t1256 = t961 * t1345;
t1255 = t963 * t1346;
t1254 = t963 * t1344;
t1202 = t1030 * t1323;
t1049 = pkin(3) ^ 2;
t1050 = pkin(2) ^ 2;
t1067 = 0.2e1 * t1018 * t1049 - t1049 + t1050 + 0.2e1 * t1291;
t992 = pkin(1) * t1357;
t951 = t1067 * t1036 + t992;
t957 = -0.2e1 * t1037 * t1268 + t1164;
t969 = t1291 + t1050 / 0.2e1 + (t1018 - 0.1e1 / 0.2e1) * t1049;
t985 = pkin(1) * t1036 - t1357;
t874 = (t991 * t1274 - t969 * t1305) * t1380 + (t1005 * t951 - t957 * t1336) * t1042 + t1002 * t1365 + t985 * t1312;
t850 = (t1371 - t874 / 0.2e1) * t1351;
t1063 = -t1005 * t1202 + t850 * t1374;
t1233 = t1036 * t1351;
t1185 = -0.2e1 * t1233;
t1112 = t930 * t1185;
t1119 = -0.2e1 * t960 * t1290;
t821 = ((pkin(2) * t1112 - t1005 * t1326) * t1042 + t850 * t1377) * t1029 + (t930 * t1119 - t1063) * t1028;
t1253 = t821 * t1351;
t875 = (t991 * t1271 + t969 * t1302) * t1380 + (t951 * t1002 + t957 * t1312) * t1042 - t1005 * t1365 + t985 * t1336;
t853 = (t1368 - t875 / 0.2e1) * t1351;
t1066 = -t1002 * t1202 + t853 * t1374;
t1111 = t933 * t1185;
t822 = ((pkin(2) * t1111 - t1002 * t1326) * t1042 + t853 * t1377) * t1029 + (t933 * t1119 - t1066) * t1028;
t1252 = t822 * t1351;
t1115 = t1290 * t1389;
t1071 = t960 * t1115;
t1124 = t1326 * t1362;
t1190 = pkin(2) * t1279;
t1205 = t1028 * t1326;
t825 = t930 * t1071 + (-t1005 * t1205 + ((-t874 + 0.2e1 * t1371) * t1029 + t930 * t1190) * t1351) * t1042 + (t1005 * t1124 + t1028 * t850) * t1377;
t1251 = t825 * t1351;
t826 = t933 * t1071 + (-t1002 * t1205 + ((-t875 + 0.2e1 * t1368) * t1029 + t933 * t1190) * t1351) * t1042 + (t1002 * t1124 + t1028 * t853) * t1377;
t1250 = t826 * t1351;
t1114 = t1289 * t1389;
t1072 = t962 * t1114;
t1123 = t1325 * t1362;
t1189 = pkin(2) * t1278;
t1204 = t1028 * t1325;
t952 = t1067 * t1038 + t992;
t958 = -0.2e1 * t1039 * t1267 + t1163;
t986 = pkin(1) * t1038 - t1357;
t876 = (t991 * t1273 - t969 * t1304) * t1379 + (t1006 * t952 - t958 * t1334) * t1044 + t1003 * t1364 + t986 * t1311;
t851 = (t1370 - t876 / 0.2e1) * t1349;
t827 = t931 * t1072 + (-t1006 * t1204 + ((-t876 + 0.2e1 * t1370) * t1029 + t931 * t1189) * t1349) * t1044 + (t1006 * t1123 + t1028 * t851) * t1376;
t1249 = t827 * t1349;
t877 = (t991 * t1270 + t969 * t1301) * t1379 + (t952 * t1003 + t958 * t1311) * t1044 - t1006 * t1364 + t986 * t1334;
t854 = (t1367 - t877 / 0.2e1) * t1349;
t828 = t934 * t1072 + (-t1003 * t1204 + ((-t877 + 0.2e1 * t1367) * t1029 + t934 * t1189) * t1349) * t1044 + (t1003 * t1123 + t1028 * t854) * t1376;
t1248 = t828 * t1349;
t1200 = t1031 * t1322;
t1062 = -t1006 * t1200 + t851 * t1373;
t1232 = t1038 * t1349;
t1184 = -0.2e1 * t1232;
t1110 = t931 * t1184;
t1118 = -0.2e1 * t962 * t1289;
t831 = ((pkin(2) * t1110 - t1006 * t1325) * t1044 + t851 * t1376) * t1029 + (t931 * t1118 - t1062) * t1028;
t1247 = t831 * t1349;
t1065 = -t1003 * t1200 + t854 * t1373;
t1109 = t934 * t1184;
t832 = ((pkin(2) * t1109 - t1003 * t1325) * t1044 + t854 * t1376) * t1029 + (t934 * t1118 - t1065) * t1028;
t1246 = t832 * t1349;
t1245 = t932 * t1341;
t1244 = t935 * t1341;
t1198 = t1032 * t1321;
t953 = t1067 * t1040 + t992;
t959 = -0.2e1 * t1041 * t1266 + t1162;
t987 = pkin(1) * t1040 - t1357;
t878 = (t991 * t1272 - t969 * t1303) * t1378 + (t1007 * t953 - t959 * t1332) * t1046 + t1004 * t1363 + t987 * t1310;
t852 = (t1369 - t878 / 0.2e1) * t1341;
t1061 = -t1007 * t1198 + t852 * t1372;
t1218 = t1040 * t1341;
t1182 = -0.2e1 * t1218;
t1108 = t932 * t1182;
t1116 = -0.2e1 * t964 * t1275;
t823 = ((pkin(2) * t1108 - t1007 * t1324) * t1046 + t852 * t1375) * t1029 + (t932 * t1116 - t1061) * t1028;
t1243 = t823 * t1341;
t879 = (t991 * t1269 + t969 * t1300) * t1378 + (t953 * t1004 + t959 * t1310) * t1046 - t1007 * t1363 + t987 * t1332;
t855 = (t1366 - t879 / 0.2e1) * t1341;
t1064 = -t1004 * t1198 + t855 * t1372;
t1107 = t935 * t1182;
t824 = ((pkin(2) * t1107 - t1004 * t1324) * t1046 + t855 * t1375) * t1029 + (t935 * t1116 - t1064) * t1028;
t1242 = t824 * t1341;
t1113 = t1275 * t1389;
t1070 = t964 * t1113;
t1122 = t1324 * t1362;
t1188 = pkin(2) * t1277;
t1203 = t1028 * t1324;
t829 = t932 * t1070 + (-t1007 * t1203 + ((-t878 + 0.2e1 * t1369) * t1029 + t932 * t1188) * t1341) * t1046 + (t1007 * t1122 + t1028 * t852) * t1375;
t1241 = t829 * t1341;
t830 = t935 * t1070 + (-t1004 * t1203 + ((-t879 + 0.2e1 * t1366) * t1029 + t935 * t1188) * t1341) * t1046 + (t1004 * t1122 + t1028 * t855) * t1375;
t1240 = t830 * t1341;
t956 = t1041 * t1017 + (pkin(1) + t984) * t1047;
t1239 = t956 * t1339;
t1238 = t956 * t1338;
t1237 = t965 * t1339;
t1236 = t965 * t1338;
t1019 = t1036 ^ 2;
t1235 = t1019 * t1350;
t1020 = t1038 ^ 2;
t1234 = t1020 * t1348;
t1231 = t960 * t1320;
t1230 = t976 * t1320;
t1229 = t960 * t1317;
t1228 = t970 * t1317;
t1227 = t973 * t1317;
t1226 = t962 * t1316;
t1225 = t978 * t1316;
t1224 = t962 * t1313;
t1223 = t1069 * t1313;
t1222 = t974 * t1313;
t1221 = -pkin(1) - t1008;
t1220 = -pkin(1) - t1009;
t1219 = -pkin(1) - t1010;
t1217 = t964 * t1308;
t1216 = t980 * t1308;
t1021 = t1040 ^ 2;
t1215 = t1021 * t1340;
t1214 = t964 * t1306;
t1213 = t1068 * t1306;
t1212 = t975 * t1306;
t1211 = t1002 * t1323;
t1210 = t1003 * t1322;
t1209 = t1004 * t1321;
t1208 = t1005 * t1323;
t1207 = t1006 * t1322;
t1206 = t1007 * t1321;
t1201 = t1030 * t1318;
t1199 = t1031 * t1314;
t1197 = t1036 * t1319;
t1196 = t1036 * t1318;
t1195 = t1038 * t1315;
t1194 = t1038 * t1314;
t1193 = t1032 * t1307;
t1192 = t1040 * t1307;
t1191 = t1001 * t1040 * t1046;
t1187 = t1351 * t1388;
t1186 = t1349 * t1387;
t1183 = t1341 * t1386;
t1181 = t1022 * t1050 + ((t1381 + qJ(3,3)) * qJ(3,3)) + t1292;
t1180 = t1024 * t1050 + ((t1381 + qJ(3,2)) * qJ(3,2)) + t1292;
t1179 = t1026 * t1050 + ((t1381 + qJ(3,1)) * qJ(3,1)) + t1292;
t1178 = t874 * t1257;
t1177 = t874 * t1256;
t1176 = t875 * t1257;
t1175 = t875 * t1256;
t1174 = t876 * t1255;
t1173 = t876 * t1254;
t1172 = t877 * t1255;
t1171 = t877 * t1254;
t1170 = t1350 * t1385;
t1169 = t1348 * t1384;
t1168 = t960 * t1261;
t1167 = t960 * t1260;
t1166 = t962 * t1259;
t1165 = t962 * t1258;
t1161 = t964 * t1239;
t1160 = t964 * t1238;
t1159 = t976 * t1233;
t1158 = t978 * t1232;
t1157 = t960 * t1230;
t1156 = t962 * t1225;
t1155 = t878 * t1237;
t1154 = t878 * t1236;
t1153 = t879 * t1237;
t1152 = t879 * t1236;
t1151 = t1340 * t1383;
t1150 = t980 * t1218;
t1149 = t964 * t1216;
t1148 = t1019 * t1229;
t1147 = t1020 * t1224;
t1146 = t1030 * t1233;
t1145 = t1030 * t1231;
t1144 = t1031 * t1232;
t1143 = t1031 * t1226;
t1142 = t961 * t1197;
t1141 = t1036 * t1229;
t1140 = t976 * t1196;
t1139 = t963 * t1195;
t1138 = t1038 * t1224;
t1137 = t978 * t1194;
t1136 = t1042 * t1229;
t1135 = t1043 * t1230;
t1134 = t1044 * t1224;
t1133 = t1045 * t1225;
t1132 = t1032 * t1218;
t1131 = t1032 * t1217;
t1130 = t980 * t1192;
t1129 = t1047 * t1216;
t1128 = t1021 * t1214;
t1127 = t965 * t1191;
t1126 = t1040 * t1214;
t1125 = t1046 * t1214;
t1121 = t1231 * t1390;
t1120 = t1226 * t1390;
t1117 = t1217 * t1390;
t1106 = t1229 * t1388;
t1105 = t1224 * t1387;
t1104 = t1214 * t1386;
t1103 = t1002 * t1159;
t1102 = t1002 * t1157;
t1101 = t1003 * t1158;
t1100 = t1003 * t1156;
t1099 = t1005 * t1159;
t1098 = t1005 * t1157;
t1097 = t1006 * t1158;
t1096 = t1006 * t1156;
t1095 = t930 * t1146;
t1094 = t931 * t1144;
t1093 = t933 * t1146;
t1092 = t934 * t1144;
t1091 = t1004 * t1150;
t1090 = t1004 * t1149;
t1089 = t1007 * t1150;
t1088 = t1007 * t1149;
t1087 = t1028 * t1145;
t1086 = t1028 * t1143;
t1085 = t1029 * t1145;
t1084 = t1029 * t1143;
t1083 = t1030 * t1140;
t1082 = t1031 * t1137;
t1081 = t1036 * t1136;
t1080 = t1038 * t1134;
t1079 = t932 * t1132;
t1078 = t935 * t1132;
t1077 = t1028 * t1131;
t1076 = t1029 * t1131;
t1075 = t1040 * t1125;
t1074 = t1004 * t1130;
t1073 = t1007 * t1130;
t1060 = t976 * (t1002 * t930 + t1005 * t933) * t1351;
t1059 = t978 * (t1003 * t931 + t1006 * t934) * t1349;
t1058 = t980 * (t1004 * t932 + t1007 * t935) * t1341;
t919 = t1002 * t1140 + t1003 * t1137 + t1074;
t920 = t1002 * t1135 + t1003 * t1133 + t1004 * t1129;
t921 = t1005 * t1140 + t1006 * t1137 + t1073;
t922 = t1005 * t1135 + t1006 * t1133 + t1007 * t1129;
t1057 = -t932 * t1091 - t931 * t1101 - t930 * t1103;
t1056 = t935 * t1091 + t934 * t1101 + t933 * t1103;
t1055 = t935 * t1090 + t934 * t1100 + t933 * t1102;
t1054 = t932 * t1089 + t931 * t1097 + t930 * t1099;
t1053 = -t935 * t1089 - t934 * t1097 - t933 * t1099;
t1052 = t932 * t1088 + t931 * t1096 + t930 * t1098;
t981 = 0.1e1 / t984 ^ 2;
t979 = 0.1e1 / t983 ^ 2;
t977 = 0.1e1 / t982 ^ 2;
t947 = t975 * t1193;
t946 = t1068 * t1193;
t945 = t974 * t1199;
t944 = t1069 * t1199;
t943 = t973 * t1201;
t942 = t970 * t1201;
t941 = (t1359 - t956 / 0.2e1) * t1000;
t940 = (t1360 - t955 / 0.2e1) * t998;
t939 = (t1361 - t954 / 0.2e1) * t996;
t938 = (t1219 * t1047 + t956) * t1000;
t937 = (t1220 * t1045 + t955) * t998;
t936 = (t1221 * t1043 + t954) * t996;
t929 = t935 ^ 2;
t928 = t934 ^ 2;
t927 = t933 ^ 2;
t926 = t932 ^ 2;
t925 = t931 ^ 2;
t924 = t930 ^ 2;
t923 = t1002 * t977 * t1005 + t1003 * t979 * t1006 + t1004 * t981 * t1007;
t918 = -pkin(2) * t1206 + t932 * t1183;
t917 = -pkin(2) * t1207 + t931 * t1186;
t916 = -pkin(2) * t1208 + t930 * t1187;
t915 = -pkin(2) * t1209 + t935 * t1183;
t914 = -pkin(2) * t1210 + t934 * t1186;
t913 = -pkin(2) * t1211 + t933 * t1187;
t912 = -t1078 + t1283;
t911 = t1078 - 0.2e1 * t1283;
t910 = -t1092 + t1284;
t909 = t1092 - 0.2e1 * t1284;
t908 = -t1093 + t1285;
t907 = t1093 - 0.2e1 * t1285;
t906 = -t1079 + t1280;
t905 = t1079 - 0.2e1 * t1280;
t904 = -t1094 + t1281;
t903 = t1094 - 0.2e1 * t1281;
t902 = -t1095 + t1282;
t901 = t1095 - 0.2e1 * t1282;
t900 = pkin(1) * t1108 - t1007 * t1286;
t899 = pkin(1) * t1107 - t1004 * t1286;
t898 = pkin(1) * t1110 - t1006 * t1287;
t897 = pkin(1) * t1109 - t1003 * t1287;
t896 = pkin(1) * t1112 - t1005 * t1288;
t895 = pkin(1) * t1111 - t1002 * t1288;
t894 = -pkin(5) * t1206 + t932 * t1117;
t893 = -pkin(5) * t1209 + t935 * t1117;
t892 = -pkin(5) * t1207 + t931 * t1120;
t891 = -pkin(5) * t1210 + t934 * t1120;
t890 = -pkin(5) * t1208 + t930 * t1121;
t889 = -pkin(5) * t1211 + t933 * t1121;
t888 = 0.2e1 * t941 * t1010 + (-pkin(1) * t956 + t1179 * t1047) * t1000;
t887 = 0.2e1 * t940 * t1009 + (-pkin(1) * t955 + t1180 * t1045) * t998;
t886 = 0.2e1 * t939 * t1008 + (-pkin(1) * t954 + t1181 * t1043) * t996;
t885 = t1047 * t1113 + ((-t956 + 0.2e1 * t1359) * t1029 + t1047 * t1188) * t1308 + t941 * t1277;
t884 = t1045 * t1114 + ((-t955 + 0.2e1 * t1360) * t1029 + t1045 * t1189) * t1316 + t940 * t1278;
t883 = t1043 * t1115 + ((-t954 + 0.2e1 * t1361) * t1029 + t1043 * t1190) * t1320 + t939 * t1279;
t882 = -0.2e1 * t1028 * t1045 * t1289 + 0.2e1 * (-t1028 * t940 - t1194 * t1358) * t1044 - 0.2e1 * t940 * t1294;
t881 = -0.2e1 * t1028 * t1047 * t1275 + 0.2e1 * (-t1028 * t941 - t1192 * t1358) * t1046 - 0.2e1 * t941 * t1293;
t880 = -0.2e1 * t1028 * t1043 * t1290 + 0.2e1 * (-t1028 * t939 - t1196 * t1358) * t1042 - 0.2e1 * t939 * t1295;
t873 = -t911 * t1029 - t935 * t1077;
t872 = t911 * t1028 - t935 * t1076;
t871 = -t909 * t1029 - t934 * t1086;
t870 = t909 * t1028 - t934 * t1084;
t869 = -t907 * t1029 - t933 * t1087;
t868 = t907 * t1028 - t933 * t1085;
t867 = -t905 * t1029 - t932 * t1077;
t866 = t905 * t1028 - t932 * t1076;
t865 = -t903 * t1029 - t931 * t1086;
t864 = t903 * t1028 - t931 * t1084;
t863 = -t901 * t1029 - t930 * t1087;
t862 = t901 * t1028 - t930 * t1085;
t861 = t935 * t1214 + t934 * t1224 + t933 * t1229;
t860 = t932 * t1214 + t931 * t1224 + t930 * t1229;
t859 = t935 * t1128 + t934 * t1147 + t933 * t1148;
t858 = t932 * t1128 + t931 * t1147 + t930 * t1148;
t857 = 0.2e1 * t935 * t1075 + 0.2e1 * t934 * t1080 + 0.2e1 * t933 * t1081;
t856 = 0.2e1 * t932 * t1075 + 0.2e1 * t931 * t1080 + 0.2e1 * t930 * t1081;
t849 = (t1219 * t935 + t879) * t1341;
t848 = (t1220 * t934 + t877) * t1349;
t847 = (t1221 * t933 + t875) * t1351;
t846 = (t1219 * t932 + t878) * t1341;
t845 = (t1220 * t931 + t876) * t1349;
t844 = (t1221 * t930 + t874) * t1351;
t843 = t1151 + t1169 + t1170;
t842 = t1019 * t1170 + t1020 * t1169 + t1021 * t1151;
t841 = 0.2e1 * t1127 * t1383 + 0.2e1 * t1139 * t1384 + 0.2e1 * t1142 * t1385;
t840 = t1042 * t1060 + t1044 * t1059 + t1046 * t1058;
t839 = t1036 * t1060 + t1038 * t1059 + t1040 * t1058;
t838 = t1064 * pkin(2) + (-pkin(1) * t879 + t1179 * t935) * t1341;
t837 = t1061 * pkin(2) + (-pkin(1) * t878 + t1179 * t932) * t1341;
t836 = t1065 * pkin(2) + (-pkin(1) * t877 + t1180 * t934) * t1349;
t835 = t1062 * pkin(2) + (-pkin(1) * t876 + t1180 * t931) * t1349;
t834 = t1066 * pkin(2) + (-pkin(1) * t875 + t1181 * t933) * t1351;
t833 = t1063 * pkin(2) + (-pkin(1) * t874 + t1181 * t930) * t1351;
t1 = [t929 * t1340 + t928 * t1348 + t927 * t1350, 0, 0, t929 * t1215 + t928 * t1234 + t927 * t1235, 0.2e1 * t929 * t1127 + 0.2e1 * t928 * t1139 + 0.2e1 * t927 * t1142, 0.2e1 * t1056, 0.2e1 * t1055, t1002 ^ 2 * t977 + t1003 ^ 2 * t979 + t1004 ^ 2 * t981, -t1056 * pkin(5) + t893 * t1244 + t891 * t1262 + t889 * t1263, -t1055 * pkin(5) + t899 * t1244 + t897 * t1262 + t895 * t1263, t869 * t1337 + t871 * t1335 + t873 * t1333 + (t1153 + t1240) * t935 + (t1172 + t1248) * t934 + (t1176 + t1250) * t933, t868 * t1337 + t870 * t1335 + t872 * t1333 + (t1152 + t1242) * t935 + (t1171 + t1246) * t934 + (t1175 + t1252) * t933, -pkin(2) * t1056 + t1244 * t915 + t1262 * t914 + t1263 * t913, (t836 * t934 + t848 * t877) * t1349 + (t834 * t933 + t847 * t875) * t1351 + (t838 * t935 + t849 * t879) * t1341 + (t912 * t1333 + t910 * t1335 + t908 * t1337) * pkin(2), 1; t843, 0, 0, t842, t841, t839, t840, t923, t1057 * pkin(5) + t894 * t1244 + t892 * t1262 + t890 * t1263, t900 * t1244 + t896 * t1263 + t898 * t1262 + (-t932 * t1090 - t931 * t1100 - t930 * t1102) * pkin(5), t932 * t1153 + t931 * t1172 + t930 * t1176 + t935 * t1241 + t934 * t1249 + t933 * t1251 + t867 * t1333 + t865 * t1335 + t863 * t1337, t932 * t1152 + t931 * t1171 + t930 * t1175 + t935 * t1243 + t934 * t1247 + t933 * t1253 + t866 * t1333 + t864 * t1335 + t862 * t1337, pkin(2) * t1057 + t1244 * t918 + t1262 * t917 + t1263 * t916, (t835 * t934 + t845 * t877) * t1349 + (t833 * t933 + t844 * t875) * t1351 + (t837 * t935 + t846 * t879) * t1341 + (t906 * t1333 + t904 * t1335 + t902 * t1337) * pkin(2), 0; t861, 0, 0, t859, t857, t919, t920, 0, -t919 * pkin(5) + (t935 * t1125 + t934 * t1134 + t933 * t1136) * t1390, -t920 * pkin(5) + (-t935 * t1126 - t934 * t1138 - t933 * t1141) * t1390, -t943 * t1337 - t945 * t1335 - t947 * t1333 + (t879 * t1213 + t885 * t1342) * t964 + (t877 * t1223 + t884 * t1352) * t962 + (t875 * t1228 + t883 * t1353) * t960, t942 * t1337 + t944 * t1335 + t946 * t1333 + (t879 * t1212 + t881 * t1342) * t964 + (t877 * t1222 + t882 * t1352) * t962 + (t875 * t1227 + t880 * t1353) * t960, -pkin(2) * t919 + t1104 * t935 + t1105 * t934 + t1106 * t933, (t877 * t937 + t887 * t934) * t1349 + (t875 * t936 + t886 * t933) * t1351 + (t879 * t938 + t888 * t935) * t1341 + (-t1002 * t1083 - t1003 * t1082 - t1032 * t1074) * pkin(2), 0; t843, 0, 0, t842, t841, t839, t840, t923, t1053 * pkin(5) + t893 * t1245 + t891 * t1264 + t889 * t1265, t899 * t1245 + t895 * t1265 + t897 * t1264 + (-t935 * t1088 - t934 * t1096 - t933 * t1098) * pkin(5), t935 * t1155 + t934 * t1174 + t933 * t1178 + t932 * t1240 + t931 * t1248 + t930 * t1250 + t873 * t1329 + t871 * t1330 + t869 * t1331, t935 * t1154 + t934 * t1173 + t933 * t1177 + t932 * t1242 + t931 * t1246 + t930 * t1252 + t872 * t1329 + t870 * t1330 + t868 * t1331, pkin(2) * t1053 + t1245 * t915 + t1264 * t914 + t1265 * t913, (t836 * t931 + t848 * t876) * t1349 + (t834 * t930 + t847 * t874) * t1351 + (t838 * t932 + t849 * t878) * t1341 + (t912 * t1329 + t910 * t1330 + t908 * t1331) * pkin(2), 0; t926 * t1340 + t925 * t1348 + t924 * t1350, 0, 0, t926 * t1215 + t925 * t1234 + t924 * t1235, 0.2e1 * t926 * t1127 + 0.2e1 * t925 * t1139 + 0.2e1 * t924 * t1142, 0.2e1 * t1054, 0.2e1 * t1052, t1005 ^ 2 * t977 + t1006 ^ 2 * t979 + t1007 ^ 2 * t981, -t1054 * pkin(5) + t894 * t1245 + t892 * t1264 + t890 * t1265, -t1052 * pkin(5) + t900 * t1245 + t898 * t1264 + t896 * t1265, t863 * t1331 + t865 * t1330 + t867 * t1329 + (t1155 + t1241) * t932 + (t1174 + t1249) * t931 + (t1178 + t1251) * t930, t862 * t1331 + t864 * t1330 + t866 * t1329 + (t1154 + t1243) * t932 + (t1173 + t1247) * t931 + (t1177 + t1253) * t930, -pkin(2) * t1054 + t1245 * t918 + t1264 * t917 + t1265 * t916, (t835 * t931 + t845 * t876) * t1349 + (t833 * t930 + t844 * t874) * t1351 + (t837 * t932 + t846 * t878) * t1341 + (t906 * t1329 + t904 * t1330 + t902 * t1331) * pkin(2), 1; t860, 0, 0, t858, t856, t921, t922, 0, -t921 * pkin(5) + (t932 * t1125 + t931 * t1134 + t930 * t1136) * t1390, -t922 * pkin(5) + (-t932 * t1126 - t931 * t1138 - t930 * t1141) * t1390, -t943 * t1331 - t945 * t1330 - t947 * t1329 + (t878 * t1213 + t885 * t1343) * t964 + (t876 * t1223 + t884 * t1354) * t962 + (t874 * t1228 + t883 * t1355) * t960, t942 * t1331 + t944 * t1330 + t946 * t1329 + (t878 * t1212 + t881 * t1343) * t964 + (t876 * t1222 + t882 * t1354) * t962 + (t1227 * t874 + t880 * t1355) * t960, -pkin(2) * t921 + t1104 * t932 + t1105 * t931 + t1106 * t930, (t876 * t937 + t887 * t931) * t1349 + (t874 * t936 + t886 * t930) * t1351 + (t878 * t938 + t888 * t932) * t1341 + (-t1005 * t1083 - t1006 * t1082 - t1032 * t1073) * pkin(2), 0; t861, 0, 0, t859, t857, t919, t920, 0, t893 * t1307 + t891 * t1314 + t889 * t1318, t899 * t1307 + t897 * t1314 + t895 * t1318, t935 * t1161 + t934 * t1166 + t933 * t1168 + t830 * t1307 + t828 * t1314 + t826 * t1318, t1160 * t935 + t1165 * t934 + t1167 * t933 + t1307 * t824 + t1314 * t832 + t1318 * t822, t1307 * t915 + t1314 * t914 + t1318 * t913, (t1045 * t836 + t848 * t955) * t998 + (t1043 * t834 + t847 * t954) * t996 + (t1047 * t838 + t849 * t956) * t1000, 0; t860, 0, 0, t858, t856, t921, t922, 0, t894 * t1307 + t892 * t1314 + t890 * t1318, t900 * t1307 + t898 * t1314 + t896 * t1318, t932 * t1161 + t931 * t1166 + t930 * t1168 + t829 * t1307 + t827 * t1314 + t825 * t1318, t1160 * t932 + t1165 * t931 + t1167 * t930 + t1307 * t823 + t1314 * t831 + t1318 * t821, t1307 * t918 + t1314 * t917 + t1318 * t916, (t1045 * t835 + t845 * t955) * t998 + (t1043 * t833 + t844 * t954) * t996 + (t1047 * t837 + t846 * t956) * t1000, 0; t1299 + t1327 + t1328, 0, 0, t1019 * t1328 + t1020 * t1327 + t1021 * t1299, 0.2e1 * t1023 * t1197 + 0.2e1 * t1025 * t1195 + 0.2e1 * t1027 * t1191, 0, 0, 0, (t1023 * t1319 + t1025 * t1315 + t1046 * t1299) * t1390, (-t1036 * t1328 - t1038 * t1327 - t1040 * t1299) * t1390, (t1000 * t885 + t1239) * t1047 + (t884 * t998 + t1259) * t1045 + (t883 * t996 + t1261) * t1043, (t1000 * t881 + t1238) * t1047 + (t882 * t998 + t1258) * t1045 + (t880 * t996 + t1260) * t1043, 0.2e1 * t1030 * t1328 + 0.2e1 * t1031 * t1327 + 0.2e1 * t1032 * t1299, (t1045 * t887 + t937 * t955) * t998 + (t1043 * t886 + t936 * t954) * t996 + (t1047 * t888 + t938 * t956) * t1000, 1;];
tau_reg  = t1;
