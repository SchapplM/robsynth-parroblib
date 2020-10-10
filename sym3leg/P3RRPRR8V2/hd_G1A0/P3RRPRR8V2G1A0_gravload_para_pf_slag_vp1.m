% Calculate Gravitation load for parallel robot
% P3RRPRR8V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:02:38
% EndTime: 2020-08-06 21:02:38
% DurationCPUTime: 0.54s
% Computational Cost: add. (666->133), mult. (1083->204), div. (27->9), fcn. (627->35), ass. (0->96)
t1206 = pkin(5) + qJ(3,1);
t1205 = pkin(5) + qJ(3,2);
t1204 = pkin(5) + qJ(3,3);
t1154 = qJ(2,3) + pkin(7);
t1203 = pkin(3) * cos(t1154);
t1155 = qJ(2,2) + pkin(7);
t1202 = pkin(3) * cos(t1155);
t1156 = qJ(2,1) + pkin(7);
t1201 = pkin(3) * cos(t1156);
t1157 = sin(pkin(7));
t1200 = pkin(3) * t1157;
t1198 = cos(pkin(7));
t1102 = -m(2) * rSges(2,1) + (-rSges(3,1) * t1198 + rSges(3,2) * t1157 - pkin(2)) * m(3);
t1199 = g(3) * t1102;
t1197 = 0.2e1 * pkin(2) * pkin(3);
t1196 = 2 * pkin(1);
t1158 = legFrame(3,3);
t1136 = sin(t1158);
t1139 = cos(t1158);
t1117 = -t1136 * g(1) + t1139 * g(2);
t1120 = t1139 * g(1) + t1136 * g(2);
t1151 = -pkin(6) - t1204;
t1142 = 0.1e1 / t1151;
t1162 = sin(qJ(1,3));
t1168 = cos(qJ(1,3));
t1195 = (-t1117 * t1168 + t1120 * t1162) * t1142;
t1159 = legFrame(2,3);
t1137 = sin(t1159);
t1140 = cos(t1159);
t1118 = -t1137 * g(1) + t1140 * g(2);
t1121 = t1140 * g(1) + t1137 * g(2);
t1152 = -pkin(6) - t1205;
t1143 = 0.1e1 / t1152;
t1164 = sin(qJ(1,2));
t1170 = cos(qJ(1,2));
t1194 = (-t1118 * t1170 + t1121 * t1164) * t1143;
t1160 = legFrame(1,3);
t1138 = sin(t1160);
t1141 = cos(t1160);
t1119 = -t1138 * g(1) + t1141 * g(2);
t1122 = t1141 * g(1) + t1138 * g(2);
t1153 = -pkin(6) - t1206;
t1144 = 0.1e1 / t1153;
t1166 = sin(qJ(1,1));
t1172 = cos(qJ(1,1));
t1193 = (-t1119 * t1172 + t1122 * t1166) * t1144;
t1116 = m(2) * rSges(2,2) + (rSges(3,1) * t1157 + t1198 * rSges(3,2)) * m(3);
t1161 = sin(qJ(2,3));
t1167 = cos(qJ(2,3));
t1189 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1182 = -t1102 * t1167 - t1116 * t1161 + t1189;
t1126 = (-rSges(2,3) - pkin(5)) * m(2) + m(1) * rSges(1,2);
t1188 = -(rSges(3,3) + t1204) * m(3) + t1126;
t1192 = t1142 * ((t1162 * t1182 + t1188 * t1168) * t1120 + (t1162 * t1188 - t1182 * t1168) * t1117);
t1163 = sin(qJ(2,2));
t1169 = cos(qJ(2,2));
t1181 = -t1102 * t1169 - t1116 * t1163 + t1189;
t1187 = -(rSges(3,3) + t1205) * m(3) + t1126;
t1191 = t1143 * ((t1164 * t1181 + t1187 * t1170) * t1121 + (t1164 * t1187 - t1181 * t1170) * t1118);
t1165 = sin(qJ(2,1));
t1171 = cos(qJ(2,1));
t1180 = -t1102 * t1171 - t1116 * t1165 + t1189;
t1186 = -(rSges(3,3) + t1206) * m(3) + t1126;
t1190 = t1144 * ((t1166 * t1180 + t1186 * t1172) * t1122 + (t1166 * t1186 - t1180 * t1172) * t1119);
t1185 = t1162 * t1117 + t1120 * t1168;
t1184 = t1164 * t1118 + t1121 * t1170;
t1183 = t1166 * t1119 + t1122 * t1172;
t1179 = pkin(2) ^ 2;
t1178 = pkin(3) ^ 2;
t1177 = 0.2e1 * qJ(2,1);
t1176 = 0.2e1 * qJ(2,2);
t1175 = 0.2e1 * qJ(2,3);
t1147 = t1171 * pkin(2);
t1146 = t1169 * pkin(2);
t1145 = t1167 * pkin(2);
t1131 = t1147 + pkin(1);
t1130 = t1146 + pkin(1);
t1129 = t1145 + pkin(1);
t1127 = t1198 * pkin(3) + pkin(2);
t1125 = 0.1e1 / (t1147 + t1201);
t1124 = 0.1e1 / (t1146 + t1202);
t1123 = 0.1e1 / (t1145 + t1203);
t1115 = g(3) * t1116;
t1114 = t1138 * t1172 + t1141 * t1166;
t1113 = t1137 * t1170 + t1140 * t1164;
t1112 = t1136 * t1168 + t1139 * t1162;
t1111 = -t1138 * t1166 + t1141 * t1172;
t1110 = -t1137 * t1164 + t1140 * t1170;
t1109 = -t1136 * t1162 + t1139 * t1168;
t1108 = t1131 * t1172 - t1166 * t1153;
t1107 = t1130 * t1170 - t1164 * t1152;
t1106 = t1129 * t1168 - t1162 * t1151;
t1105 = t1166 * t1131 + t1172 * t1153;
t1104 = t1164 * t1130 + t1170 * t1152;
t1103 = t1162 * t1129 + t1168 * t1151;
t1 = [-t1109 * t1192 - t1110 * t1191 - t1111 * t1190 - m(4) * g(1) + ((-t1105 * t1138 + t1108 * t1141 + t1111 * t1201) * t1193 + (-t1104 * t1137 + t1107 * t1140 + t1110 * t1202) * t1194 + (-t1103 * t1136 + t1106 * t1139 + t1109 * t1203) * t1195) * m(3); -t1112 * t1192 - t1113 * t1191 - t1114 * t1190 - m(4) * g(2) + ((t1105 * t1141 + t1108 * t1138 + t1114 * t1201) * t1193 + (t1104 * t1140 + t1107 * t1137 + t1113 * t1202) * t1194 + (t1103 * t1139 + t1106 * t1136 + t1112 * t1203) * t1195) * m(3); -(t1165 * t1127 + t1171 * t1200) / (t1127 * t1171 - t1165 * t1200) * t1190 + t1125 * ((-t1183 * t1102 + t1115) * t1165 + (t1183 * t1116 + t1199) * t1171) - (t1163 * t1127 + t1169 * t1200) / (t1127 * t1169 - t1163 * t1200) * t1191 + t1124 * ((-t1184 * t1102 + t1115) * t1163 + (t1184 * t1116 + t1199) * t1169) - (t1161 * t1127 + t1167 * t1200) / (t1127 * t1167 - t1161 * t1200) * t1192 + t1123 * ((-t1185 * t1102 + t1115) * t1161 + (t1185 * t1116 + t1199) * t1167) - m(4) * g(3) + ((sin(t1177 + pkin(7)) * t1197 + t1179 * sin(t1177) + t1178 * sin(0.2e1 * t1156) + (sin(t1156) * pkin(3) + pkin(2) * t1165) * t1196) * t1125 * t1193 + (sin(t1176 + pkin(7)) * t1197 + t1179 * sin(t1176) + t1178 * sin(0.2e1 * t1155) + (sin(t1155) * pkin(3) + pkin(2) * t1163) * t1196) * t1124 * t1194 + (sin(t1175 + pkin(7)) * t1197 + t1179 * sin(t1175) + t1178 * sin(0.2e1 * t1154) + (sin(t1154) * pkin(3) + pkin(2) * t1161) * t1196) * t1123 * t1195) * m(3) / 0.2e1;];
taugX  = t1;
