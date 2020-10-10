% Calculate Gravitation load for parallel robot
% P3RPRRR9V1G3A0
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
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:56:42
% EndTime: 2020-08-06 18:56:43
% DurationCPUTime: 0.75s
% Computational Cost: add. (654->176), mult. (891->281), div. (42->10), fcn. (636->26), ass. (0->129)
t1178 = 2 * pkin(3);
t1098 = cos(pkin(7));
t1177 = 0.2e1 * t1098 ^ 2;
t1117 = m(2) + m(3);
t1101 = pkin(5) + qJ(2,1);
t1100 = pkin(5) + qJ(2,2);
t1099 = pkin(5) + qJ(2,3);
t1111 = cos(qJ(3,3));
t1094 = t1111 ^ 2;
t1176 = t1094 * pkin(3);
t1113 = cos(qJ(3,2));
t1095 = t1113 ^ 2;
t1175 = t1095 * pkin(3);
t1115 = cos(qJ(3,1));
t1096 = t1115 ^ 2;
t1174 = t1096 * pkin(3);
t1173 = t1111 * pkin(2);
t1172 = t1113 * pkin(2);
t1171 = t1115 * pkin(2);
t1102 = legFrame(3,2);
t1078 = sin(t1102);
t1081 = cos(t1102);
t1041 = t1078 * g(1) + t1081 * g(2);
t1091 = pkin(7) + qJ(3,3);
t1071 = cos(t1091);
t1065 = 0.1e1 / t1071;
t1068 = sin(t1091);
t1044 = t1081 * g(1) - t1078 * g(2);
t1106 = sin(qJ(1,3));
t1112 = cos(qJ(1,3));
t1132 = -g(3) * t1106 + t1044 * t1112;
t1170 = ((-mrSges(3,1) * t1041 + t1132 * mrSges(3,2)) * t1071 + t1068 * (t1132 * mrSges(3,1) + mrSges(3,2) * t1041)) * t1065;
t1103 = legFrame(2,2);
t1079 = sin(t1103);
t1082 = cos(t1103);
t1042 = t1079 * g(1) + t1082 * g(2);
t1092 = pkin(7) + qJ(3,2);
t1072 = cos(t1092);
t1066 = 0.1e1 / t1072;
t1069 = sin(t1092);
t1045 = t1082 * g(1) - t1079 * g(2);
t1108 = sin(qJ(1,2));
t1114 = cos(qJ(1,2));
t1131 = -g(3) * t1108 + t1045 * t1114;
t1169 = ((-mrSges(3,1) * t1042 + t1131 * mrSges(3,2)) * t1072 + t1069 * (t1131 * mrSges(3,1) + mrSges(3,2) * t1042)) * t1066;
t1104 = legFrame(1,2);
t1080 = sin(t1104);
t1083 = cos(t1104);
t1043 = t1080 * g(1) + t1083 * g(2);
t1093 = pkin(7) + qJ(3,1);
t1073 = cos(t1093);
t1067 = 0.1e1 / t1073;
t1070 = sin(t1093);
t1046 = t1083 * g(1) - t1080 * g(2);
t1110 = sin(qJ(1,1));
t1116 = cos(qJ(1,1));
t1130 = -g(3) * t1110 + t1046 * t1116;
t1168 = ((-mrSges(3,1) * t1043 + t1130 * mrSges(3,2)) * t1073 + t1070 * (t1130 * mrSges(3,1) + mrSges(3,2) * t1043)) * t1067;
t1087 = pkin(6) + t1099;
t1075 = 0.1e1 / t1087;
t1167 = (-t1112 * g(3) - t1044 * t1106) * t1075;
t1088 = pkin(6) + t1100;
t1076 = 0.1e1 / t1088;
t1166 = (-t1114 * g(3) - t1045 * t1108) * t1076;
t1089 = pkin(6) + t1101;
t1077 = 0.1e1 / t1089;
t1165 = (-t1116 * g(3) - t1046 * t1110) * t1077;
t1105 = sin(qJ(3,3));
t1119 = pkin(2) / 0.2e1;
t1164 = (t1111 * pkin(3) + t1119) * t1105;
t1107 = sin(qJ(3,2));
t1163 = (t1113 * pkin(3) + t1119) * t1107;
t1109 = sin(qJ(3,1));
t1162 = (t1115 * pkin(3) + t1119) * t1109;
t1060 = pkin(1) * t1117 + mrSges(1,1);
t1056 = g(3) * t1060;
t1143 = mrSges(2,3) + mrSges(3,3) - mrSges(1,2);
t1124 = m(2) * qJ(2,3) + t1099 * m(3) + t1143;
t1097 = sin(pkin(7));
t1142 = -m(3) * pkin(2) - mrSges(2,1);
t1129 = -(-mrSges(3,1) * t1111 + mrSges(3,2) * t1105 + t1142) * t1098 - (t1105 * mrSges(3,1) + mrSges(3,2) * t1111 + mrSges(2,2)) * t1097;
t1161 = t1075 * (t1056 * t1112 + (t1124 * t1106 + t1129 * t1112) * g(3) + (-t1124 * t1112 - (-t1060 - t1129) * t1106) * t1044);
t1125 = m(2) * qJ(2,2) + t1100 * m(3) + t1143;
t1128 = -(-mrSges(3,1) * t1113 + mrSges(3,2) * t1107 + t1142) * t1098 - (t1107 * mrSges(3,1) + mrSges(3,2) * t1113 + mrSges(2,2)) * t1097;
t1160 = t1076 * (t1056 * t1114 + (t1125 * t1108 + t1128 * t1114) * g(3) + (-t1125 * t1114 - (-t1060 - t1128) * t1108) * t1045);
t1126 = m(2) * qJ(2,1) + t1101 * m(3) + t1143;
t1127 = -(-mrSges(3,1) * t1115 + mrSges(3,2) * t1109 + t1142) * t1098 - (t1109 * mrSges(3,1) + mrSges(3,2) * t1115 + mrSges(2,2)) * t1097;
t1159 = t1077 * (t1056 * t1116 + (t1126 * t1110 + t1127 * t1116) * g(3) + (-t1126 * t1116 - (-t1060 - t1127) * t1110) * t1046);
t1158 = t1078 * t1112;
t1157 = t1079 * t1114;
t1156 = t1080 * t1116;
t1155 = t1081 * t1112;
t1154 = t1082 * t1114;
t1153 = t1083 * t1116;
t1152 = t1097 * t1105;
t1151 = t1097 * t1107;
t1150 = t1097 * t1109;
t1074 = pkin(1) * t1097;
t1149 = t1111 * (-t1105 * pkin(3) + t1074);
t1148 = t1113 * (-t1107 * pkin(3) + t1074);
t1147 = t1115 * (-t1109 * pkin(3) + t1074);
t1146 = pkin(1) * t1112 + t1106 * t1087;
t1145 = pkin(1) * t1114 + t1108 * t1088;
t1144 = pkin(1) * t1116 + t1110 * t1089;
t1141 = t1065 * t1161;
t1140 = t1066 * t1160;
t1139 = t1067 * t1159;
t1138 = 0.1e1 / (t1098 * t1111 - t1152) * t1167;
t1137 = 0.1e1 / (t1098 * t1113 - t1151) * t1166;
t1136 = 0.1e1 / (t1098 * t1115 - t1150) * t1165;
t1135 = t1112 * t1152;
t1134 = t1114 * t1151;
t1133 = t1116 * t1150;
t1123 = pkin(2) * t1135 + (t1135 * t1178 - t1146) * t1111;
t1122 = pkin(2) * t1134 + (t1134 * t1178 - t1145) * t1113;
t1121 = pkin(2) * t1133 + (t1133 * t1178 - t1144) * t1115;
t1120 = 1 / pkin(3);
t1118 = -pkin(3) / 0.2e1;
t1064 = t1098 * pkin(2) + pkin(1);
t1049 = t1174 + t1171 / 0.2e1 + t1118;
t1048 = t1175 + t1172 / 0.2e1 + t1118;
t1047 = t1176 + t1173 / 0.2e1 + t1118;
t1028 = pkin(1) * t1109 + (-pkin(3) + t1171 + 0.2e1 * t1174) * t1097;
t1027 = pkin(1) * t1107 + (-pkin(3) + t1172 + 0.2e1 * t1175) * t1097;
t1026 = pkin(1) * t1105 + (-pkin(3) + t1173 + 0.2e1 * t1176) * t1097;
t1025 = t1144 * t1150 + (t1096 - 0.1e1) * t1116 * pkin(3);
t1024 = t1145 * t1151 + (t1095 - 0.1e1) * t1114 * pkin(3);
t1023 = t1146 * t1152 + (t1094 - 0.1e1) * t1112 * pkin(3);
t1 = [(t1080 * t1070 + t1073 * t1153) * t1139 + (t1079 * t1069 + t1072 * t1154) * t1140 + (t1078 * t1068 + t1071 * t1155) * t1141 - g(1) * m(4) + (t1078 * t1170 + t1079 * t1169 + t1080 * t1168) * t1120 + (((t1049 * t1153 + t1080 * t1162) * t1177 + (t1080 * t1028 - t1121 * t1083) * t1098 - t1025 * t1083 + t1080 * t1147) * t1136 + ((t1048 * t1154 + t1079 * t1163) * t1177 + (t1079 * t1027 - t1122 * t1082) * t1098 - t1024 * t1082 + t1079 * t1148) * t1137 + ((t1047 * t1155 + t1078 * t1164) * t1177 + (t1078 * t1026 - t1123 * t1081) * t1098 - t1023 * t1081 + t1078 * t1149) * t1138) * t1117; (t1083 * t1070 - t1073 * t1156) * t1139 + (t1082 * t1069 - t1072 * t1157) * t1140 + (t1081 * t1068 - t1071 * t1158) * t1141 - g(2) * m(4) + (t1081 * t1170 + t1082 * t1169 + t1083 * t1168) * t1120 + (((-t1049 * t1156 + t1083 * t1162) * t1177 + (t1083 * t1028 + t1121 * t1080) * t1098 + t1025 * t1080 + t1083 * t1147) * t1136 + ((-t1048 * t1157 + t1082 * t1163) * t1177 + (t1082 * t1027 + t1122 * t1079) * t1098 + t1024 * t1079 + t1082 * t1148) * t1137 + ((-t1047 * t1158 + t1081 * t1164) * t1177 + (t1081 * t1026 + t1123 * t1078) * t1098 + t1023 * t1078 + t1081 * t1149) * t1138) * t1117; -t1106 * t1161 - t1108 * t1160 - t1110 * t1159 - g(3) * m(4) + ((t1089 * t1116 + (-pkin(3) * t1073 - t1064) * t1110) * t1165 + (t1088 * t1114 + (-pkin(3) * t1072 - t1064) * t1108) * t1166 + (t1087 * t1112 + (-pkin(3) * t1071 - t1064) * t1106) * t1167) * t1117;];
taugX  = t1;
