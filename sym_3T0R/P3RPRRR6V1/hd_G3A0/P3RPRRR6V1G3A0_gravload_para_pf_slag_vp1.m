% Calculate Gravitation load for parallel robot
% P3RPRRR6V1G3A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:41:27
% EndTime: 2020-08-06 18:41:29
% DurationCPUTime: 1.46s
% Computational Cost: add. (762->200), mult. (936->281), div. (51->10), fcn. (765->53), ass. (0->121)
t1219 = sin(pkin(7));
t1237 = -pkin(6) - pkin(5);
t1179 = t1237 * t1219 - pkin(1);
t1231 = cos(qJ(1,3));
t1170 = t1179 * t1231;
t1220 = cos(pkin(7));
t1225 = sin(qJ(1,3));
t1275 = t1225 * t1237;
t1276 = t1225 * t1219;
t1301 = t1170 + pkin(2) * t1276 - (pkin(2) * t1231 - t1275) * t1220;
t1233 = cos(qJ(1,2));
t1171 = t1179 * t1233;
t1227 = sin(qJ(1,2));
t1271 = t1227 * t1237;
t1272 = t1227 * t1219;
t1300 = t1171 + pkin(2) * t1272 - (pkin(2) * t1233 - t1271) * t1220;
t1235 = cos(qJ(1,1));
t1172 = t1179 * t1235;
t1229 = sin(qJ(1,1));
t1267 = t1229 * t1237;
t1268 = t1229 * t1219;
t1299 = t1172 + pkin(2) * t1268 - (pkin(2) * t1235 - t1267) * t1220;
t1223 = legFrame(1,2);
t1200 = sin(t1223);
t1203 = cos(t1223);
t1169 = t1203 * g(1) - t1200 * g(2);
t1234 = cos(qJ(3,1));
t1190 = t1234 * pkin(3) + pkin(2);
t1197 = t1220 * pkin(1);
t1175 = 0.1e1 / (t1197 + t1190);
t1212 = qJ(1,1) + pkin(7);
t1193 = sin(t1212);
t1196 = cos(t1212);
t1209 = t1235 * pkin(1);
t1228 = sin(qJ(3,1));
t1246 = t1234 * rSges(3,1) - t1228 * rSges(3,2);
t1242 = pkin(2) + t1246;
t1292 = rSges(3,3) + pkin(5);
t1260 = g(3) * t1292;
t1288 = t1229 * pkin(1);
t1283 = t1175 * (-m(1) * (t1169 * (-t1229 * rSges(1,1) - t1235 * rSges(1,2)) - g(3) * (t1235 * rSges(1,1) - t1229 * rSges(1,2))) - m(2) * (t1169 * (-t1193 * rSges(2,1) - t1196 * rSges(2,2) - t1288) - g(3) * (t1196 * rSges(2,1) - t1193 * rSges(2,2) + t1209)) - m(3) * (-t1169 * t1288 - g(3) * t1209 + (-g(3) * t1242 + t1169 * t1292) * t1196 + (-t1169 * t1242 - t1260) * t1193));
t1298 = t1283 / 0.2e1;
t1222 = legFrame(2,2);
t1199 = sin(t1222);
t1202 = cos(t1222);
t1168 = t1202 * g(1) - t1199 * g(2);
t1232 = cos(qJ(3,2));
t1189 = t1232 * pkin(3) + pkin(2);
t1174 = 0.1e1 / (t1197 + t1189);
t1211 = qJ(1,2) + pkin(7);
t1192 = sin(t1211);
t1195 = cos(t1211);
t1207 = t1233 * pkin(1);
t1226 = sin(qJ(3,2));
t1248 = t1232 * rSges(3,1) - t1226 * rSges(3,2);
t1243 = pkin(2) + t1248;
t1289 = t1227 * pkin(1);
t1285 = t1174 * (-m(1) * (t1168 * (-t1227 * rSges(1,1) - t1233 * rSges(1,2)) - g(3) * (t1233 * rSges(1,1) - t1227 * rSges(1,2))) - m(2) * (t1168 * (-t1192 * rSges(2,1) - t1195 * rSges(2,2) - t1289) - g(3) * (t1195 * rSges(2,1) - t1192 * rSges(2,2) + t1207)) - m(3) * (-t1168 * t1289 - g(3) * t1207 + (-g(3) * t1243 + t1168 * t1292) * t1195 + (-t1168 * t1243 - t1260) * t1192));
t1297 = t1285 / 0.2e1;
t1221 = legFrame(3,2);
t1198 = sin(t1221);
t1201 = cos(t1221);
t1167 = t1201 * g(1) - t1198 * g(2);
t1230 = cos(qJ(3,3));
t1188 = t1230 * pkin(3) + pkin(2);
t1173 = 0.1e1 / (t1197 + t1188);
t1210 = qJ(1,3) + pkin(7);
t1191 = sin(t1210);
t1194 = cos(t1210);
t1205 = t1231 * pkin(1);
t1224 = sin(qJ(3,3));
t1250 = t1230 * rSges(3,1) - t1224 * rSges(3,2);
t1244 = pkin(2) + t1250;
t1290 = t1225 * pkin(1);
t1287 = t1173 * (-m(1) * (t1167 * (-t1225 * rSges(1,1) - t1231 * rSges(1,2)) - g(3) * (t1231 * rSges(1,1) - t1225 * rSges(1,2))) - m(2) * (t1167 * (-t1191 * rSges(2,1) - t1194 * rSges(2,2) - t1290) - g(3) * (t1194 * rSges(2,1) - t1191 * rSges(2,2) + t1205)) - m(3) * (-t1167 * t1290 - g(3) * t1205 + (-g(3) * t1244 + t1167 * t1292) * t1194 + (-t1167 * t1244 - t1260) * t1191));
t1296 = t1287 / 0.2e1;
t1295 = 0.2e1 * pkin(2);
t1294 = 0.2e1 * t1237;
t1293 = -m(2) - m(3);
t1291 = m(3) / pkin(3);
t1286 = t1173 / t1224;
t1284 = t1174 / t1226;
t1282 = t1175 / t1228;
t1281 = t1188 * t1231;
t1280 = t1189 * t1233;
t1279 = t1190 * t1235;
t1278 = t1224 * t1198;
t1277 = t1224 * t1201;
t1274 = t1226 * t1199;
t1273 = t1226 * t1202;
t1270 = t1228 * t1200;
t1269 = t1228 * t1203;
t1266 = -pkin(7) + qJ(3,1);
t1265 = -pkin(7) + qJ(3,2);
t1264 = -pkin(7) + qJ(3,3);
t1263 = pkin(7) + qJ(3,1);
t1262 = pkin(7) + qJ(3,2);
t1261 = pkin(7) + qJ(3,3);
t1259 = pkin(3) * (-t1231 * t1220 + t1276) * t1230 ^ 2;
t1258 = pkin(3) * (-t1233 * t1220 + t1272) * t1232 ^ 2;
t1257 = pkin(3) * (-t1235 * t1220 + t1268) * t1234 ^ 2;
t1164 = t1198 * g(1) + t1201 * g(2);
t1256 = t1293 * t1164 * t1286;
t1165 = t1199 * g(1) + t1202 * g(2);
t1255 = t1293 * t1165 * t1284;
t1166 = t1200 * g(1) + t1203 * g(2);
t1254 = t1293 * t1166 * t1282;
t1152 = t1164 * t1250 + (-g(3) * t1191 + t1167 * t1194) * (-rSges(3,1) * t1224 - rSges(3,2) * t1230);
t1253 = t1152 * ((-t1275 + t1281) * t1220 - t1170 - t1188 * t1276) * t1286;
t1153 = t1165 * t1248 + (-g(3) * t1192 + t1168 * t1195) * (-rSges(3,1) * t1226 - rSges(3,2) * t1232);
t1252 = t1153 * ((-t1271 + t1280) * t1220 - t1171 - t1189 * t1272) * t1284;
t1154 = t1166 * t1246 + (-g(3) * t1193 + t1169 * t1196) * (-rSges(3,1) * t1228 - rSges(3,2) * t1234);
t1251 = t1154 * ((-t1267 + t1279) * t1220 - t1172 - t1190 * t1268) * t1282;
t1187 = t1197 + pkin(2);
t1186 = -t1223 + t1212;
t1185 = t1223 + t1212;
t1184 = -t1222 + t1211;
t1183 = t1222 + t1211;
t1182 = -t1221 + t1210;
t1181 = t1221 + t1210;
t1 = [(-t1203 * t1257 + (pkin(3) * t1270 - t1299 * t1203) * t1234 + t1187 * t1270) * t1254 + (-t1202 * t1258 + (pkin(3) * t1274 - t1300 * t1202) * t1232 + t1187 * t1274) * t1255 + (-t1201 * t1259 + (pkin(3) * t1278 - t1301 * t1201) * t1230 + t1187 * t1278) * t1256 - m(4) * g(1) + (-sin(t1185) - sin(t1186)) * t1298 + (-sin(t1183) - sin(t1184)) * t1297 + (-sin(t1181) - sin(t1182)) * t1296 + (t1201 * t1253 + t1202 * t1252 + t1203 * t1251) * t1291; (t1200 * t1257 + (pkin(3) * t1269 + t1299 * t1200) * t1234 + t1187 * t1269) * t1254 + (t1199 * t1258 + (pkin(3) * t1273 + t1300 * t1199) * t1232 + t1187 * t1273) * t1255 + (t1198 * t1259 + (pkin(3) * t1277 + t1301 * t1198) * t1230 + t1187 * t1277) * t1256 - m(4) * g(2) + (cos(t1186) - cos(t1185)) * t1298 + (cos(t1184) - cos(t1183)) * t1297 + (cos(t1182) - cos(t1181)) * t1296 + (-t1198 * t1253 - t1199 * t1252 - t1200 * t1251) * t1291; -t1196 * t1283 - t1234 * ((t1190 * t1229 + t1235 * t1237) * t1220 - t1179 * t1229 + t1219 * t1279) * t1254 - t1195 * t1285 - t1232 * ((t1189 * t1227 + t1233 * t1237) * t1220 - t1179 * t1227 + t1219 * t1280) * t1255 - t1194 * t1287 - t1230 * ((t1188 * t1225 + t1231 * t1237) * t1220 - t1179 * t1225 + t1219 * t1281) * t1256 - m(4) * g(3) + (-(t1196 * t1294 + 0.2e1 * t1288 + t1193 * t1295 + (sin(qJ(1,1) - t1266) + sin(qJ(1,1) + t1263)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t1228 * t1295 + (sin(t1263) + sin(t1266)) * pkin(1)) * t1154 - (t1195 * t1294 + 0.2e1 * t1289 + t1192 * t1295 + (sin(qJ(1,2) - t1265) + sin(qJ(1,2) + t1262)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t1226 * t1295 + (sin(t1262) + sin(t1265)) * pkin(1)) * t1153 - (t1194 * t1294 + 0.2e1 * t1290 + t1191 * t1295 + (sin(qJ(1,3) - t1264) + sin(qJ(1,3) + t1261)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t1224 * t1295 + (sin(t1261) + sin(t1264)) * pkin(1)) * t1152) * t1291;];
taugX  = t1;
