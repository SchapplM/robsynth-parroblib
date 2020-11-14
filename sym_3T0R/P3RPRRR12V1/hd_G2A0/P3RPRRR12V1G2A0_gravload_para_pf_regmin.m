% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR12V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x14]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR12V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:55
% EndTime: 2020-08-06 18:24:55
% DurationCPUTime: 0.66s
% Computational Cost: add. (405->123), mult. (747->229), div. (102->7), fcn. (735->18), ass. (0->111)
t1132 = sin(qJ(1,1));
t1138 = cos(qJ(1,1));
t1126 = legFrame(1,2);
t1116 = sin(t1126);
t1119 = cos(t1126);
t1149 = t1119 * g(1) - t1116 * g(2);
t1146 = g(3) * t1138 + t1149 * t1132;
t1131 = sin(qJ(3,1));
t1113 = t1131 * pkin(3) + qJ(2,1);
t1110 = 0.1e1 / t1113;
t1170 = t1138 * t1110;
t1164 = t1146 * t1170;
t1130 = sin(qJ(1,2));
t1136 = cos(qJ(1,2));
t1125 = legFrame(2,2);
t1115 = sin(t1125);
t1118 = cos(t1125);
t1150 = t1118 * g(1) - t1115 * g(2);
t1147 = g(3) * t1136 + t1150 * t1130;
t1129 = sin(qJ(3,2));
t1112 = t1129 * pkin(3) + qJ(2,2);
t1109 = 0.1e1 / t1112;
t1171 = t1136 * t1109;
t1166 = t1147 * t1171;
t1128 = sin(qJ(1,3));
t1134 = cos(qJ(1,3));
t1124 = legFrame(3,2);
t1114 = sin(t1124);
t1117 = cos(t1124);
t1151 = t1117 * g(1) - t1114 * g(2);
t1148 = g(3) * t1134 + t1151 * t1128;
t1127 = sin(qJ(3,3));
t1111 = t1127 * pkin(3) + qJ(2,3);
t1108 = 0.1e1 / t1111;
t1172 = t1134 * t1108;
t1168 = t1148 * t1172;
t1205 = g(3) * t1128 - t1151 * t1134;
t1204 = g(3) * t1130 - t1150 * t1136;
t1203 = g(3) * t1132 - t1149 * t1138;
t1202 = pkin(3) * t1134;
t1201 = pkin(3) * t1136;
t1200 = pkin(3) * t1138;
t1196 = qJ(2,2) * t1136;
t1195 = qJ(2,3) * t1134;
t1194 = t1138 * qJ(2,1);
t1121 = 0.1e1 / t1127;
t1193 = t1205 * t1121;
t1122 = 0.1e1 / t1129;
t1192 = t1204 * t1122;
t1123 = 0.1e1 / t1131;
t1191 = t1203 * t1123;
t1190 = t1108 * t1128;
t1189 = t1109 * t1130;
t1188 = t1110 * t1132;
t1187 = t1114 * t1121;
t1133 = cos(qJ(3,3));
t1186 = t1114 * t1133;
t1185 = t1115 * t1122;
t1135 = cos(qJ(3,2));
t1184 = t1115 * t1135;
t1183 = t1116 * t1123;
t1137 = cos(qJ(3,1));
t1182 = t1116 * t1137;
t1181 = t1117 * t1121;
t1180 = t1117 * t1133;
t1179 = t1118 * t1122;
t1178 = t1118 * t1135;
t1177 = t1119 * t1123;
t1176 = t1119 * t1137;
t1084 = -g(3) * (-t1128 * pkin(1) + t1195) - t1151 * (t1134 * pkin(1) + t1128 * qJ(2,3));
t1175 = t1128 * t1084;
t1085 = -g(3) * (-t1130 * pkin(1) + t1196) - t1150 * (t1136 * pkin(1) + t1130 * qJ(2,2));
t1174 = t1130 * t1085;
t1086 = -g(3) * (-t1132 * pkin(1) + t1194) - t1149 * (t1138 * pkin(1) + t1132 * qJ(2,1));
t1173 = t1132 * t1086;
t1169 = t1148 * t1190;
t1167 = t1147 * t1189;
t1165 = t1146 * t1188;
t1163 = t1114 * t1190;
t1162 = t1117 * t1190;
t1161 = t1115 * t1189;
t1160 = t1118 * t1189;
t1159 = t1116 * t1188;
t1158 = t1119 * t1188;
t1157 = t1127 * t1169;
t1156 = t1133 * t1169;
t1155 = t1129 * t1167;
t1154 = t1135 * t1167;
t1153 = t1131 * t1165;
t1152 = t1137 * t1165;
t1145 = t1170 * t1203 + t1171 * t1204 + t1172 * t1205;
t1144 = t1166 + t1168 + t1164;
t1143 = t1159 * t1203 + t1161 * t1204 + t1163 * t1205;
t1142 = t1158 * t1203 + t1160 * t1204 + t1162 * t1205;
t1141 = t1146 * t1159 + t1147 * t1161 + t1148 * t1163;
t1140 = t1146 * t1158 + t1147 * t1160 + t1148 * t1162;
t1139 = 0.1e1 / pkin(3);
t1120 = pkin(1) + pkin(5) + pkin(6);
t1107 = -t1120 * t1130 + t1196;
t1106 = -t1120 * t1128 + t1195;
t1105 = -t1120 * t1132 + t1194;
t1104 = t1116 * g(1) + t1119 * g(2);
t1103 = t1115 * g(1) + t1118 * g(2);
t1102 = t1114 * g(1) + t1117 * g(2);
t1083 = t1104 * t1131 - t1137 * t1203;
t1082 = t1104 * t1137 + t1131 * t1203;
t1081 = t1103 * t1129 - t1135 * t1204;
t1080 = t1103 * t1135 + t1129 * t1204;
t1079 = t1102 * t1127 - t1133 * t1205;
t1078 = t1102 * t1133 + t1127 * t1205;
t1 = [0, t1142, t1140, -t1142, -t1140, (t1119 * t1173 - ((pkin(3) * t1182 - t1105 * t1119) * t1131 + (t1137 - 0.1e1) * (t1137 + 0.1e1) * t1119 * t1200 + qJ(2,1) * t1182) * t1191) * t1110 + (t1118 * t1174 - ((pkin(3) * t1184 - t1107 * t1118) * t1129 + (t1135 - 0.1e1) * (t1135 + 0.1e1) * t1118 * t1201 + qJ(2,2) * t1184) * t1192) * t1109 + (t1117 * t1175 - ((pkin(3) * t1186 - t1106 * t1117) * t1127 + (t1133 - 0.1e1) * (t1133 + 0.1e1) * t1117 * t1202 + qJ(2,3) * t1186) * t1193) * t1108, 0, 0, 0, 0, 0, -t1117 * t1157 - t1118 * t1155 - t1119 * t1153 + (-t1079 * t1187 - t1081 * t1185 - t1083 * t1183) * t1139, -t1117 * t1156 - t1118 * t1154 - t1119 * t1152 + (-t1078 * t1187 - t1080 * t1185 - t1082 * t1183) * t1139, -g(1); 0, -t1143, -t1141, t1143, t1141, (-t1116 * t1173 - ((pkin(3) * t1176 + t1105 * t1116) * t1131 + (-t1137 ^ 2 + 0.1e1) * t1116 * t1200 + qJ(2,1) * t1176) * t1191) * t1110 + (-t1115 * t1174 - ((pkin(3) * t1178 + t1107 * t1115) * t1129 + (-t1135 ^ 2 + 0.1e1) * t1115 * t1201 + qJ(2,2) * t1178) * t1192) * t1109 + (-t1114 * t1175 - ((pkin(3) * t1180 + t1106 * t1114) * t1127 + (-t1133 ^ 2 + 0.1e1) * t1114 * t1202 + qJ(2,3) * t1180) * t1193) * t1108, 0, 0, 0, 0, 0, t1114 * t1157 + t1115 * t1155 + t1116 * t1153 + (-t1079 * t1181 - t1081 * t1179 - t1083 * t1177) * t1139, t1114 * t1156 + t1115 * t1154 + t1116 * t1152 + (-t1078 * t1181 - t1080 * t1179 - t1082 * t1177) * t1139, -g(2); 0, t1145, t1144, -t1145, -t1144, (t1138 * t1086 - (t1132 * t1113 + t1120 * t1138) * t1203) * t1110 + (t1136 * t1085 - (t1130 * t1112 + t1120 * t1136) * t1204) * t1109 + (t1134 * t1084 - (t1128 * t1111 + t1120 * t1134) * t1205) * t1108, 0, 0, 0, 0, 0, -t1127 * t1168 - t1129 * t1166 - t1131 * t1164, -t1133 * t1168 - t1135 * t1166 - t1137 * t1164, -g(3);];
tau_reg  = t1;
