% Calculate vector of centrifugal and coriolis load on the joints for
% P4RRRRR2G1A0
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
% koppelP [4x3]
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
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:24:39
% EndTime: 2020-08-07 17:24:46
% DurationCPUTime: 7.18s
% Computational Cost: add. (30388->470), mult. (26404->873), div. (12868->23), fcn. (23096->58), ass. (0->422)
t1346 = 2 * pkin(1);
t1116 = sin(qJ(3,4));
t1118 = cos(qJ(3,4));
t1117 = sin(qJ(2,4));
t1335 = pkin(1) * t1117;
t1004 = t1118 * (-mrSges(3,2) * t1335 + Ifges(3,6)) - t1116 * (mrSges(3,1) * t1335 - Ifges(3,5));
t1093 = 0.1e1 / t1118 ^ 2;
t1092 = 0.1e1 / t1118;
t1135 = xDP(3);
t1151 = 0.1e1 / pkin(2);
t1269 = t1135 * t1151;
t1224 = t1092 * t1269;
t1090 = 0.1e1 / t1117;
t1154 = 1 / pkin(1);
t1299 = t1090 * t1154;
t1225 = t1116 * t1299;
t1297 = t1093 * t1151;
t1119 = cos(qJ(2,4));
t1334 = pkin(1) * t1119;
t1237 = (pkin(2) * t1118 + t1334) * t1297;
t1170 = t1225 * t1237;
t1166 = t1135 * t1170;
t1139 = xP(4);
t1088 = sin(t1139);
t1089 = cos(t1139);
t1142 = koppelP(4,2);
t1146 = koppelP(4,1);
t1024 = -t1088 * t1142 + t1089 * t1146;
t1134 = xDP(4);
t1136 = xDP(2);
t1008 = t1024 * t1134 + t1136;
t1032 = 0.1e1 / (sin(qJ(2,4) + qJ(3,4)) + sin(qJ(2,4) - qJ(3,4)));
t1268 = t1151 * t1154;
t1241 = t1032 * t1268;
t1084 = qJ(1,4) + legFrame(4,3);
t1076 = qJ(2,4) + t1084;
t1055 = qJ(3,4) + t1076;
t1056 = -qJ(3,4) + t1076;
t1337 = -2 * pkin(1);
t996 = sin(t1084) * t1337 + (-sin(t1056) - sin(t1055)) * pkin(2);
t1203 = t996 * t1241;
t968 = t1008 * t1203;
t1020 = t1088 * t1146 + t1089 * t1142;
t1137 = xDP(1);
t1012 = -t1020 * t1134 + t1137;
t997 = cos(t1084) * t1337 + (-cos(t1055) - cos(t1056)) * pkin(2);
t1202 = t997 * t1241;
t969 = t1012 * t1202;
t928 = t968 + t969 - t1166;
t1195 = t1092 * t1225;
t1054 = cos(t1076);
t1232 = t1054 * t1299;
t1053 = sin(t1076);
t1233 = t1053 * t1299;
t961 = t1008 * t1233 + t1012 * t1232 + t1135 * t1195;
t911 = t928 + t961;
t1179 = -0.2e1 * t911 * t1224;
t1204 = 0.4e1 * Ifges(3,4) * t1269;
t1111 = t1135 ^ 2;
t1286 = t1111 / pkin(2) ^ 2;
t1217 = t1116 * t1286;
t1254 = Ifges(3,5) * t1286;
t1120 = Ifges(3,2) - Ifges(3,1);
t1284 = t1116 * t1120;
t1318 = t1118 * t911;
t1162 = t1204 * t1318 + Ifges(3,4) * t1179 + (t1093 * t1254 + t1179 * t1284) * t1118 - Ifges(3,6) * t1093 * t1217;
t1194 = t1119 * t1224;
t1174 = t911 * t1194;
t1091 = t1118 ^ 2;
t1155 = t1118 * t1091;
t1094 = 0.1e1 / t1155;
t1193 = t1094 * t1217;
t1121 = mrSges(2,2) - mrSges(3,3);
t1283 = t1119 * t1121;
t1209 = t1286 / 0.2e1;
t1319 = t1117 * (t1093 * t1209 + (t961 + t928 / 0.2e1) * t928);
t1153 = pkin(1) ^ 2;
t1285 = t1116 * t1117;
t1213 = t1135 * t1285;
t1327 = pkin(2) * t1091;
t1258 = t1119 * t1327;
t1150 = pkin(2) ^ 2;
t902 = t1150 * t911 * t1155;
t909 = t969 / 0.2e1 + t968 / 0.2e1 - t1166 / 0.2e1 + t961;
t876 = (((-pkin(1) * t911 * t1285 + t1092 * t1135) * t1118 + pkin(1) * t1194) * t1094 * t1269 + ((t902 + t909 * t1258 * t1346 + (-pkin(1) * t1092 * t1213 + t1153 * t961) * t1118) * t961 + (t902 + (t911 * t1258 - t1213) * pkin(1)) * t928) * t1297) * t1299;
t1287 = t1111 * t1151;
t884 = ((-t1118 * t961 * t1334 - t911 * t1327) * t1092 * t961 - pkin(2) * t928 * t1318 - t1094 * t1287) * t1299;
t1208 = mrSges(3,2) * t1116 - mrSges(2,1);
t1033 = t1208 * t1334;
t1072 = t1121 * t1335;
t1080 = mrSges(3,1) * t1334;
t1336 = (Ifges(3,1) + Ifges(2,3));
t1175 = (m(2) + m(3)) * t1153 + Ifges(1,3) + t1336;
t1282 = t1120 * t1091;
t1323 = Ifges(3,4) * t1116;
t980 = t1282 + 0.2e1 * (t1080 + t1323) * t1118 - 0.2e1 * t1033 - 0.2e1 * t1072 + t1175;
t1183 = t1282 + t1336;
t1266 = 0.2e1 * t1323;
t984 = (t1080 + t1266) * t1118 - t1033 - t1072 + t1183;
t1345 = t1004 * t1193 - t876 * t984 - t884 * t980 + t1162 + ((-mrSges(3,1) * t1319 - mrSges(3,2) * t1174) * t1118 + (-mrSges(3,1) * t1174 + mrSges(3,2) * t1319) * t1116 + (-mrSges(2,1) * t1117 - t1283) * t928 * t909) * t1346;
t1016 = t1118 * t1266 + t1183;
t1044 = Ifges(3,5) * t1116 + Ifges(3,6) * t1118;
t960 = t961 ^ 2;
t1344 = -t1016 * t876 + t1044 * t1193 - t884 * t984 + (t1283 + (mrSges(3,1) * t1118 - t1208) * t1117) * t960 * pkin(1) + t1162;
t1122 = sin(qJ(3,3));
t1128 = cos(qJ(3,3));
t1123 = sin(qJ(2,3));
t1333 = pkin(1) * t1123;
t1005 = t1128 * (-mrSges(3,2) * t1333 + Ifges(3,6)) - t1122 * (mrSges(3,1) * t1333 - Ifges(3,5));
t1100 = 0.1e1 / t1128 ^ 2;
t1099 = 0.1e1 / t1128;
t1220 = t1099 * t1269;
t1095 = 0.1e1 / t1123;
t1296 = t1095 * t1154;
t1223 = t1122 * t1296;
t1292 = t1100 * t1151;
t1129 = cos(qJ(2,3));
t1330 = pkin(1) * t1129;
t1236 = (pkin(2) * t1128 + t1330) * t1292;
t1169 = t1223 * t1236;
t1165 = t1135 * t1169;
t1143 = koppelP(3,2);
t1147 = koppelP(3,1);
t1025 = -t1088 * t1143 + t1089 * t1147;
t1009 = t1025 * t1134 + t1136;
t1035 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t1240 = t1035 * t1268;
t1085 = qJ(1,3) + legFrame(3,3);
t1077 = qJ(2,3) + t1085;
t1066 = qJ(3,3) + t1077;
t1067 = -qJ(3,3) + t1077;
t998 = sin(t1085) * t1337 + (-sin(t1067) - sin(t1066)) * pkin(2);
t1201 = t998 * t1240;
t970 = t1009 * t1201;
t1021 = t1088 * t1147 + t1089 * t1143;
t1013 = -t1021 * t1134 + t1137;
t1001 = cos(t1085) * t1337 + (-cos(t1066) - cos(t1067)) * pkin(2);
t1198 = t1001 * t1240;
t973 = t1013 * t1198;
t930 = t970 + t973 - t1165;
t1192 = t1099 * t1223;
t1063 = cos(t1077);
t1228 = t1063 * t1296;
t1060 = sin(t1077);
t1231 = t1060 * t1296;
t965 = t1009 * t1231 + t1013 * t1228 + t1135 * t1192;
t922 = t930 + t965;
t1178 = -0.2e1 * t922 * t1220;
t1216 = t1122 * t1286;
t1278 = t1120 * t1122;
t1314 = t1128 * t922;
t1161 = t1204 * t1314 + (t1100 * t1254 + t1178 * t1278) * t1128 - Ifges(3,6) * t1100 * t1216 + Ifges(3,4) * t1178;
t1189 = t1129 * t1220;
t1173 = t922 * t1189;
t1098 = t1128 ^ 2;
t1156 = t1128 * t1098;
t1101 = 0.1e1 / t1156;
t1188 = t1101 * t1216;
t1275 = t1121 * t1129;
t1317 = t1123 * (t1100 * t1209 + (t965 + t930 / 0.2e1) * t930);
t1272 = t1122 * t1123;
t1212 = t1135 * t1272;
t1326 = pkin(2) * t1098;
t1257 = t1129 * t1326;
t913 = t1150 * t922 * t1156;
t916 = t973 / 0.2e1 + t970 / 0.2e1 - t1165 / 0.2e1 + t965;
t877 = (((-pkin(1) * t922 * t1272 + t1099 * t1135) * t1128 + pkin(1) * t1189) * t1101 * t1269 + ((t913 + t916 * t1257 * t1346 + (-pkin(1) * t1099 * t1212 + t1153 * t965) * t1128) * t965 + (t913 + (t922 * t1257 - t1212) * pkin(1)) * t930) * t1292) * t1296;
t885 = ((-t1128 * t965 * t1330 - t922 * t1326) * t1099 * t965 - pkin(2) * t930 * t1314 - t1101 * t1287) * t1296;
t1207 = mrSges(3,2) * t1122 - mrSges(2,1);
t1038 = t1207 * t1330;
t1073 = t1121 * t1333;
t1081 = mrSges(3,1) * t1330;
t1281 = t1120 * t1098;
t1322 = Ifges(3,4) * t1122;
t981 = t1281 + 0.2e1 * (t1081 + t1322) * t1128 - 0.2e1 * t1038 - 0.2e1 * t1073 + t1175;
t1182 = t1281 + t1336;
t1265 = 0.2e1 * t1322;
t985 = (t1081 + t1265) * t1128 - t1038 - t1073 + t1182;
t1343 = t1005 * t1188 - t877 * t985 - t885 * t981 + t1161 + ((-mrSges(3,1) * t1317 - mrSges(3,2) * t1173) * t1128 + (-mrSges(3,1) * t1173 + mrSges(3,2) * t1317) * t1122 + (-mrSges(2,1) * t1123 - t1275) * t930 * t916) * t1346;
t1124 = sin(qJ(3,2));
t1130 = cos(qJ(3,2));
t1125 = sin(qJ(2,2));
t1332 = pkin(1) * t1125;
t1006 = t1130 * (-mrSges(3,2) * t1332 + Ifges(3,6)) - t1124 * (mrSges(3,1) * t1332 - Ifges(3,5));
t1104 = 0.1e1 / t1130 ^ 2;
t1103 = 0.1e1 / t1130;
t1219 = t1103 * t1269;
t1096 = 0.1e1 / t1125;
t1295 = t1096 * t1154;
t1222 = t1124 * t1295;
t1290 = t1104 * t1151;
t1131 = cos(qJ(2,2));
t1329 = pkin(1) * t1131;
t1235 = (pkin(2) * t1130 + t1329) * t1290;
t1168 = t1222 * t1235;
t1164 = t1135 * t1168;
t1144 = koppelP(2,2);
t1148 = koppelP(2,1);
t1026 = -t1088 * t1144 + t1089 * t1148;
t1010 = t1026 * t1134 + t1136;
t1036 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t1239 = t1036 * t1268;
t1086 = qJ(1,2) + legFrame(2,3);
t1078 = qJ(2,2) + t1086;
t1068 = qJ(3,2) + t1078;
t1069 = -qJ(3,2) + t1078;
t999 = sin(t1086) * t1337 + (-sin(t1069) - sin(t1068)) * pkin(2);
t1200 = t999 * t1239;
t971 = t1010 * t1200;
t1022 = t1088 * t1148 + t1089 * t1144;
t1014 = -t1022 * t1134 + t1137;
t1002 = cos(t1086) * t1337 + (-cos(t1068) - cos(t1069)) * pkin(2);
t1197 = t1002 * t1239;
t974 = t1014 * t1197;
t931 = t971 + t974 - t1164;
t1191 = t1103 * t1222;
t1064 = cos(t1078);
t1227 = t1064 * t1295;
t1061 = sin(t1078);
t1230 = t1061 * t1295;
t966 = t1010 * t1230 + t1014 * t1227 + t1135 * t1191;
t923 = t931 + t966;
t1177 = -0.2e1 * t923 * t1219;
t1215 = t1124 * t1286;
t1277 = t1120 * t1124;
t1313 = t1130 * t923;
t1160 = t1204 * t1313 + (t1104 * t1254 + t1177 * t1277) * t1130 - Ifges(3,6) * t1104 * t1215 + Ifges(3,4) * t1177;
t1187 = t1131 * t1219;
t1172 = t923 * t1187;
t1102 = t1130 ^ 2;
t1157 = t1130 * t1102;
t1105 = 0.1e1 / t1157;
t1186 = t1105 * t1215;
t1274 = t1121 * t1131;
t1316 = t1125 * (t1104 * t1209 + (t966 + t931 / 0.2e1) * t931);
t1271 = t1124 * t1125;
t1211 = t1135 * t1271;
t1325 = pkin(2) * t1102;
t1256 = t1131 * t1325;
t914 = t1150 * t923 * t1157;
t917 = t974 / 0.2e1 + t971 / 0.2e1 - t1164 / 0.2e1 + t966;
t878 = (((-pkin(1) * t923 * t1271 + t1103 * t1135) * t1130 + pkin(1) * t1187) * t1105 * t1269 + ((t914 + t917 * t1256 * t1346 + (-pkin(1) * t1103 * t1211 + t1153 * t966) * t1130) * t966 + (t914 + (t923 * t1256 - t1211) * pkin(1)) * t931) * t1290) * t1295;
t886 = ((-t1130 * t966 * t1329 - t923 * t1325) * t1103 * t966 - pkin(2) * t931 * t1313 - t1105 * t1287) * t1295;
t1206 = mrSges(3,2) * t1124 - mrSges(2,1);
t1039 = t1206 * t1329;
t1074 = t1121 * t1332;
t1082 = mrSges(3,1) * t1329;
t1280 = t1120 * t1102;
t1321 = Ifges(3,4) * t1124;
t982 = t1280 + 0.2e1 * (t1082 + t1321) * t1130 - 0.2e1 * t1039 - 0.2e1 * t1074 + t1175;
t1181 = t1280 + t1336;
t1264 = 0.2e1 * t1321;
t986 = (t1082 + t1264) * t1130 - t1039 - t1074 + t1181;
t1342 = t1006 * t1186 - t878 * t986 - t886 * t982 + t1160 + ((-mrSges(3,1) * t1316 - mrSges(3,2) * t1172) * t1130 + (-mrSges(3,1) * t1172 + mrSges(3,2) * t1316) * t1124 + (-mrSges(2,1) * t1125 - t1274) * t931 * t917) * t1346;
t1126 = sin(qJ(3,1));
t1132 = cos(qJ(3,1));
t1127 = sin(qJ(2,1));
t1331 = pkin(1) * t1127;
t1007 = t1132 * (-mrSges(3,2) * t1331 + Ifges(3,6)) - t1126 * (mrSges(3,1) * t1331 - Ifges(3,5));
t1108 = 0.1e1 / t1132 ^ 2;
t1107 = 0.1e1 / t1132;
t1218 = t1107 * t1269;
t1097 = 0.1e1 / t1127;
t1294 = t1097 * t1154;
t1221 = t1126 * t1294;
t1288 = t1108 * t1151;
t1133 = cos(qJ(2,1));
t1328 = pkin(1) * t1133;
t1234 = (pkin(2) * t1132 + t1328) * t1288;
t1167 = t1221 * t1234;
t1163 = t1135 * t1167;
t1145 = koppelP(1,2);
t1149 = koppelP(1,1);
t1027 = -t1088 * t1145 + t1089 * t1149;
t1011 = t1027 * t1134 + t1136;
t1087 = qJ(1,1) + legFrame(1,3);
t1079 = qJ(2,1) + t1087;
t1070 = qJ(3,1) + t1079;
t1071 = -qJ(3,1) + t1079;
t1000 = sin(t1087) * t1337 + (-sin(t1071) - sin(t1070)) * pkin(2);
t1037 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t1238 = t1037 * t1268;
t1199 = t1000 * t1238;
t972 = t1011 * t1199;
t1023 = t1088 * t1149 + t1089 * t1145;
t1015 = -t1023 * t1134 + t1137;
t1003 = cos(t1087) * t1337 + (-cos(t1070) - cos(t1071)) * pkin(2);
t1196 = t1003 * t1238;
t975 = t1015 * t1196;
t932 = t972 + t975 - t1163;
t1190 = t1107 * t1221;
t1065 = cos(t1079);
t1226 = t1065 * t1294;
t1062 = sin(t1079);
t1229 = t1062 * t1294;
t967 = t1011 * t1229 + t1015 * t1226 + t1135 * t1190;
t924 = t932 + t967;
t1176 = -0.2e1 * t924 * t1218;
t1214 = t1126 * t1286;
t1276 = t1120 * t1126;
t1312 = t1132 * t924;
t1159 = t1204 * t1312 + (t1108 * t1254 + t1176 * t1276) * t1132 - Ifges(3,6) * t1108 * t1214 + Ifges(3,4) * t1176;
t1185 = t1133 * t1218;
t1171 = t924 * t1185;
t1106 = t1132 ^ 2;
t1158 = t1132 * t1106;
t1109 = 0.1e1 / t1158;
t1184 = t1109 * t1214;
t1273 = t1121 * t1133;
t1315 = t1127 * (t1108 * t1209 + (t967 + t932 / 0.2e1) * t932);
t1270 = t1126 * t1127;
t1210 = t1135 * t1270;
t1324 = pkin(2) * t1106;
t1255 = t1133 * t1324;
t915 = t1150 * t924 * t1158;
t918 = t975 / 0.2e1 + t972 / 0.2e1 - t1163 / 0.2e1 + t967;
t879 = (((-pkin(1) * t924 * t1270 + t1107 * t1135) * t1132 + pkin(1) * t1185) * t1109 * t1269 + ((t915 + t918 * t1255 * t1346 + (-pkin(1) * t1107 * t1210 + t1153 * t967) * t1132) * t967 + (t915 + (t924 * t1255 - t1210) * pkin(1)) * t932) * t1288) * t1294;
t887 = ((-t1132 * t967 * t1328 - t924 * t1324) * t1107 * t967 - pkin(2) * t932 * t1312 - t1109 * t1287) * t1294;
t1205 = mrSges(3,2) * t1126 - mrSges(2,1);
t1040 = t1205 * t1328;
t1075 = t1121 * t1331;
t1083 = mrSges(3,1) * t1328;
t1279 = t1120 * t1106;
t1320 = Ifges(3,4) * t1126;
t983 = t1279 + 0.2e1 * (t1083 + t1320) * t1132 - 0.2e1 * t1040 - 0.2e1 * t1075 + t1175;
t1180 = t1279 + t1336;
t1263 = 0.2e1 * t1320;
t987 = (t1083 + t1263) * t1132 - t1040 - t1075 + t1180;
t1341 = t1007 * t1184 - t879 * t987 - t887 * t983 + t1159 + ((-mrSges(3,1) * t1315 - mrSges(3,2) * t1171) * t1132 + (-mrSges(3,1) * t1171 + mrSges(3,2) * t1315) * t1126 + (-mrSges(2,1) * t1127 - t1273) * t932 * t918) * t1346;
t1017 = t1128 * t1265 + t1182;
t1046 = Ifges(3,5) * t1122 + Ifges(3,6) * t1128;
t962 = t965 ^ 2;
t1340 = -t1017 * t877 + t1046 * t1188 - t885 * t985 + (t1275 + (mrSges(3,1) * t1128 - t1207) * t1123) * t962 * pkin(1) + t1161;
t1018 = t1130 * t1264 + t1181;
t1047 = Ifges(3,5) * t1124 + Ifges(3,6) * t1130;
t963 = t966 ^ 2;
t1339 = -t1018 * t878 + t1047 * t1186 - t886 * t986 + (t1274 + (mrSges(3,1) * t1130 - t1206) * t1125) * t963 * pkin(1) + t1160;
t1019 = t1132 * t1263 + t1180;
t1048 = Ifges(3,5) * t1126 + Ifges(3,6) * t1132;
t964 = t967 ^ 2;
t1338 = -t1019 * t879 + t1048 * t1184 - t887 * t987 + (t1273 + (mrSges(3,1) * t1132 - t1205) * t1127) * t964 * pkin(1) + t1159;
t1110 = t1134 ^ 2;
t1311 = t1032 * t1151;
t1310 = t1035 * t1151;
t1309 = t1036 * t1151;
t1308 = t1037 * t1151;
t1307 = t1053 * t1090;
t1306 = t1054 * t1090;
t1305 = t1060 * t1095;
t1304 = t1061 * t1096;
t1303 = t1062 * t1097;
t1302 = t1063 * t1095;
t1301 = t1064 * t1096;
t1300 = t1065 * t1097;
t1298 = t1092 * t1151;
t1293 = t1099 * t1151;
t1291 = t1103 * t1151;
t1289 = t1107 * t1151;
t1267 = t1154 * t1110;
t1262 = t960 * t1334;
t1261 = t962 * t1330;
t1260 = t963 * t1329;
t1259 = t964 * t1328;
t1253 = t996 * t1311;
t1252 = t997 * t1311;
t1251 = t998 * t1310;
t1250 = t999 * t1309;
t1249 = t1000 * t1308;
t1248 = t1001 * t1310;
t1247 = t1002 * t1309;
t1246 = t1003 * t1308;
t1245 = t1016 * t1311;
t1244 = t1017 * t1310;
t1243 = t1018 * t1309;
t1242 = t1019 * t1308;
t1141 = mrSges(4,1);
t1140 = mrSges(4,2);
t979 = (-t1023 * t1065 + t1027 * t1062) * t1294;
t978 = (-t1022 * t1064 + t1026 * t1061) * t1295;
t977 = (-t1021 * t1063 + t1025 * t1060) * t1296;
t976 = (-t1020 * t1054 + t1024 * t1053) * t1299;
t959 = t1048 * t1289 + (-t1019 * t1234 + t1107 * t987) * t1221;
t958 = t1047 * t1291 + (-t1018 * t1235 + t1103 * t986) * t1222;
t957 = t1046 * t1293 + (-t1017 * t1236 + t1099 * t985) * t1223;
t956 = t1044 * t1298 + (-t1016 * t1237 + t1092 * t984) * t1225;
t955 = (t1003 * t1242 + t987 * t1300) * t1154;
t954 = (t1002 * t1243 + t986 * t1301) * t1154;
t953 = (t1001 * t1244 + t985 * t1302) * t1154;
t952 = (t1000 * t1242 + t987 * t1303) * t1154;
t951 = (t999 * t1243 + t986 * t1304) * t1154;
t950 = (t998 * t1244 + t985 * t1305) * t1154;
t949 = (t1000 * t1027 - t1003 * t1023) * t1238;
t948 = (-t1002 * t1022 + t1026 * t999) * t1239;
t947 = (-t1001 * t1021 + t1025 * t998) * t1240;
t946 = (t997 * t1245 + t984 * t1306) * t1154;
t945 = (t996 * t1245 + t984 * t1307) * t1154;
t944 = (-t1020 * t997 + t1024 * t996) * t1241;
t943 = (t987 * t1246 + t983 * t1300) * t1154;
t942 = (t986 * t1247 + t982 * t1301) * t1154;
t941 = (t985 * t1248 + t981 * t1302) * t1154;
t940 = (t987 * t1249 + t983 * t1303) * t1154;
t939 = (t986 * t1250 + t982 * t1304) * t1154;
t938 = (t985 * t1251 + t981 * t1305) * t1154;
t937 = t1007 * t1289 + (t1107 * t983 - t987 * t1234) * t1221;
t936 = t1006 * t1291 + (t1103 * t982 - t986 * t1235) * t1222;
t935 = t1005 * t1293 + (t1099 * t981 - t985 * t1236) * t1223;
t934 = (t984 * t1252 + t980 * t1306) * t1154;
t933 = (t984 * t1253 + t980 * t1307) * t1154;
t929 = t1004 * t1298 + (t1092 * t980 - t984 * t1237) * t1225;
t927 = t1019 * t949 + t979 * t987;
t926 = t1018 * t948 + t978 * t986;
t925 = t1017 * t947 + t977 * t985;
t921 = t924 ^ 2;
t920 = t923 ^ 2;
t919 = t922 ^ 2;
t912 = t1016 * t944 + t976 * t984;
t910 = t911 ^ 2;
t895 = t949 * t987 + t979 * t983;
t894 = t948 * t986 + t978 * t982;
t893 = t947 * t985 + t977 * t981;
t892 = t944 * t984 + t976 * t980;
t1 = [t1110 * (t1088 * t1140 - t1089 * t1141) + t1345 * t1232 + t1343 * t1228 + t1342 * t1227 + t1341 * t1226 + t1344 * t1202 + t1340 * t1198 + t1339 * t1197 + t1338 * t1196 + (-(t955 * t1246 + t943 * t1300) * t1027 - (t955 * t1249 + t943 * t1303) * t1023 - (t954 * t1247 + t942 * t1301) * t1026 - (t954 * t1250 + t942 * t1304) * t1022 - (t953 * t1248 + t941 * t1302) * t1025 - (t953 * t1251 + t941 * t1305) * t1021 - (t946 * t1252 + t934 * t1306) * t1024 - (t946 * t1253 + t934 * t1307) * t1020) * t1267; -t1110 * (t1088 * t1141 + t1089 * t1140) + t1345 * t1233 + t1343 * t1231 + t1342 * t1230 + t1341 * t1229 + t1344 * t1203 + t1340 * t1201 + t1339 * t1200 + t1338 * t1199 + (-(t952 * t1246 + t940 * t1300) * t1027 - (t952 * t1249 + t940 * t1303) * t1023 - (t951 * t1247 + t939 * t1301) * t1026 - (t951 * t1250 + t939 * t1304) * t1022 - (t950 * t1248 + t938 * t1302) * t1025 - (t950 * t1251 + t938 * t1305) * t1021 - (t945 * t1252 + t933 * t1306) * t1024 - (t945 * t1253 + t933 * t1307) * t1020) * t1267; ((mrSges(3,2) * t1262 + t910 * t1284) * t1118 + mrSges(3,1) * t1116 * t1262 + (-0.2e1 * t1091 + 0.1e1) * t910 * Ifges(3,4) + Ifges(3,3) * t1193 - t1004 * t884 - t1044 * t876) * t1298 + ((mrSges(3,2) * t1261 + t919 * t1278) * t1128 + mrSges(3,1) * t1122 * t1261 + (-0.2e1 * t1098 + 0.1e1) * t919 * Ifges(3,4) + Ifges(3,3) * t1188 - t1005 * t885 - t1046 * t877) * t1293 + ((mrSges(3,2) * t1260 + t920 * t1277) * t1130 + mrSges(3,1) * t1124 * t1260 + (-0.2e1 * t1102 + 0.1e1) * t920 * Ifges(3,4) + Ifges(3,3) * t1186 - t1006 * t886 - t1047 * t878) * t1291 + ((mrSges(3,2) * t1259 + t921 * t1276) * t1132 + mrSges(3,1) * t1126 * t1259 + (-0.2e1 * t1106 + 0.1e1) * t921 * Ifges(3,4) + Ifges(3,3) * t1184 - t1007 * t887 - t1048 * t879) * t1289 + t1345 * t1195 + t1343 * t1192 + t1342 * t1191 + t1341 * t1190 - t1344 * t1170 - t1340 * t1169 - t1339 * t1168 - t1338 * t1167 + (-(t956 * t1252 + t929 * t1306) * t1024 - (t956 * t1253 + t929 * t1307) * t1020 - (t959 * t1246 + t937 * t1300) * t1027 - (t959 * t1249 + t937 * t1303) * t1023 - (t958 * t1247 + t936 * t1301) * t1026 - (t958 * t1250 + t936 * t1304) * t1022 - (t957 * t1248 + t935 * t1302) * t1025 - (t957 * t1251 + t935 * t1305) * t1021) * t1267; t1341 * t979 + t1342 * t978 + t1343 * t977 + t1345 * t976 + t1338 * t949 + t1339 * t948 + t1340 * t947 + t1344 * t944 + (-(t927 * t1246 + t895 * t1300) * t1027 - (t927 * t1249 + t895 * t1303) * t1023 - (t926 * t1247 + t894 * t1301) * t1026 - (t926 * t1250 + t894 * t1304) * t1022 - (t925 * t1248 + t893 * t1302) * t1025 - (t925 * t1251 + t893 * t1305) * t1021 - (t912 * t1252 + t892 * t1306) * t1024 - (t912 * t1253 + t892 * t1307) * t1020) * t1267;];
taucX  = t1;
