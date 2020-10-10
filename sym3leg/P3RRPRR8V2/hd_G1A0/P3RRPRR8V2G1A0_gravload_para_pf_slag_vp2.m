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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
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

function taugX = P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:02:59
% EndTime: 2020-08-06 21:03:00
% DurationCPUTime: 0.56s
% Computational Cost: add. (630->131), mult. (795->199), div. (27->9), fcn. (555->35), ass. (0->96)
t1134 = -qJ(3,1) - pkin(5);
t1133 = -qJ(3,2) - pkin(5);
t1132 = -qJ(3,3) - pkin(5);
t1125 = qJ(2,3) + pkin(7);
t1176 = pkin(3) * cos(t1125);
t1126 = qJ(2,2) + pkin(7);
t1175 = pkin(3) * cos(t1126);
t1127 = qJ(2,1) + pkin(7);
t1174 = pkin(3) * cos(t1127);
t1128 = sin(pkin(7));
t1173 = pkin(3) * t1128;
t1172 = cos(pkin(7));
t1171 = 0.2e1 * pkin(2) * pkin(3);
t1170 = 2 * pkin(1);
t1129 = legFrame(3,3);
t1110 = sin(t1129);
t1113 = cos(t1129);
t1091 = -t1110 * g(1) + t1113 * g(2);
t1094 = t1113 * g(1) + t1110 * g(2);
t1122 = -pkin(6) + t1132;
t1116 = 0.1e1 / t1122;
t1136 = sin(qJ(1,3));
t1142 = cos(qJ(1,3));
t1169 = (-t1142 * t1091 + t1136 * t1094) * t1116;
t1130 = legFrame(2,3);
t1111 = sin(t1130);
t1114 = cos(t1130);
t1092 = -t1111 * g(1) + t1114 * g(2);
t1095 = t1114 * g(1) + t1111 * g(2);
t1123 = -pkin(6) + t1133;
t1117 = 0.1e1 / t1123;
t1138 = sin(qJ(1,2));
t1144 = cos(qJ(1,2));
t1168 = (-t1144 * t1092 + t1138 * t1095) * t1117;
t1131 = legFrame(1,3);
t1112 = sin(t1131);
t1115 = cos(t1131);
t1093 = -t1112 * g(1) + t1115 * g(2);
t1096 = t1115 * g(1) + t1112 * g(2);
t1124 = -pkin(6) + t1134;
t1118 = 0.1e1 / t1124;
t1140 = sin(qJ(1,1));
t1146 = cos(qJ(1,1));
t1167 = (-t1146 * t1093 + t1140 * t1096) * t1118;
t1089 = -m(3) * pkin(2) - mrSges(3,1) * t1172 + mrSges(3,2) * t1128 - mrSges(2,1);
t1100 = t1128 * mrSges(3,1) + t1172 * mrSges(3,2) + mrSges(2,2);
t1101 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t1135 = sin(qJ(2,3));
t1141 = cos(qJ(2,3));
t1156 = t1089 * t1141 + t1100 * t1135 - t1101;
t1163 = -m(2) * pkin(5) + mrSges(1,2) - mrSges(2,3) - mrSges(3,3);
t1159 = t1132 * m(3) + t1163;
t1166 = t1116 * ((-t1136 * t1156 + t1159 * t1142) * t1094 + (t1136 * t1159 + t1156 * t1142) * t1091);
t1137 = sin(qJ(2,2));
t1143 = cos(qJ(2,2));
t1155 = t1089 * t1143 + t1100 * t1137 - t1101;
t1158 = t1133 * m(3) + t1163;
t1165 = t1117 * ((-t1138 * t1155 + t1158 * t1144) * t1095 + (t1138 * t1158 + t1155 * t1144) * t1092);
t1139 = sin(qJ(2,1));
t1145 = cos(qJ(2,1));
t1154 = t1089 * t1145 + t1100 * t1139 - t1101;
t1157 = t1134 * m(3) + t1163;
t1164 = t1118 * ((-t1140 * t1154 + t1157 * t1146) * t1096 + (t1140 * t1157 + t1154 * t1146) * t1093);
t1162 = -t1136 * t1091 - t1094 * t1142;
t1161 = -t1138 * t1092 - t1095 * t1144;
t1160 = -t1140 * t1093 - t1096 * t1146;
t1153 = pkin(2) ^ 2;
t1152 = pkin(3) ^ 2;
t1151 = 0.2e1 * qJ(2,1);
t1150 = 0.2e1 * qJ(2,2);
t1149 = 0.2e1 * qJ(2,3);
t1121 = t1145 * pkin(2);
t1120 = t1143 * pkin(2);
t1119 = t1141 * pkin(2);
t1105 = t1121 + pkin(1);
t1104 = t1120 + pkin(1);
t1103 = t1119 + pkin(1);
t1102 = t1172 * pkin(3) + pkin(2);
t1099 = 0.1e1 / (t1121 + t1174);
t1098 = 0.1e1 / (t1120 + t1175);
t1097 = 0.1e1 / (t1119 + t1176);
t1090 = g(3) * t1100;
t1088 = g(3) * t1089;
t1087 = t1112 * t1146 + t1115 * t1140;
t1086 = t1111 * t1144 + t1114 * t1138;
t1085 = t1110 * t1142 + t1113 * t1136;
t1084 = -t1112 * t1140 + t1115 * t1146;
t1083 = -t1111 * t1138 + t1114 * t1144;
t1082 = -t1110 * t1136 + t1113 * t1142;
t1081 = t1105 * t1146 - t1140 * t1124;
t1080 = t1104 * t1144 - t1138 * t1123;
t1079 = t1103 * t1142 - t1136 * t1122;
t1078 = t1140 * t1105 + t1146 * t1124;
t1077 = t1138 * t1104 + t1144 * t1123;
t1076 = t1136 * t1103 + t1142 * t1122;
t1 = [-t1082 * t1166 - t1083 * t1165 - t1084 * t1164 - g(1) * m(4) + ((-t1078 * t1112 + t1081 * t1115 + t1084 * t1174) * t1167 + (-t1077 * t1111 + t1080 * t1114 + t1083 * t1175) * t1168 + (-t1076 * t1110 + t1079 * t1113 + t1082 * t1176) * t1169) * m(3); -t1085 * t1166 - t1086 * t1165 - t1087 * t1164 - g(2) * m(4) + ((t1078 * t1115 + t1081 * t1112 + t1087 * t1174) * t1167 + (t1077 * t1114 + t1080 * t1111 + t1086 * t1175) * t1168 + (t1076 * t1113 + t1079 * t1110 + t1085 * t1176) * t1169) * m(3); -(t1139 * t1102 + t1145 * t1173) / (t1102 * t1145 - t1139 * t1173) * t1164 + t1099 * ((t1160 * t1089 + t1090) * t1139 - (t1160 * t1100 - t1088) * t1145) - (t1137 * t1102 + t1143 * t1173) / (t1102 * t1143 - t1137 * t1173) * t1165 + t1098 * ((t1161 * t1089 + t1090) * t1137 - (t1161 * t1100 - t1088) * t1143) - (t1135 * t1102 + t1141 * t1173) / (t1102 * t1141 - t1135 * t1173) * t1166 + t1097 * ((t1162 * t1089 + t1090) * t1135 - (t1162 * t1100 - t1088) * t1141) - g(3) * m(4) + ((sin(t1151 + pkin(7)) * t1171 + t1153 * sin(t1151) + t1152 * sin(0.2e1 * t1127) + (sin(t1127) * pkin(3) + pkin(2) * t1139) * t1170) * t1099 * t1167 + (sin(t1150 + pkin(7)) * t1171 + t1153 * sin(t1150) + t1152 * sin(0.2e1 * t1126) + (sin(t1126) * pkin(3) + pkin(2) * t1137) * t1170) * t1098 * t1168 + (sin(t1149 + pkin(7)) * t1171 + t1153 * sin(t1149) + t1152 * sin(0.2e1 * t1125) + (sin(t1125) * pkin(3) + pkin(2) * t1135) * t1170) * t1097 * t1169) * m(3) / 0.2e1;];
taugX  = t1;
