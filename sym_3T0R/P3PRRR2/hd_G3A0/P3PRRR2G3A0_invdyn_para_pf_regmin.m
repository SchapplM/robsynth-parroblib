% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRR2G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x8]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRR2G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:05
% EndTime: 2020-03-09 21:20:06
% DurationCPUTime: 0.93s
% Computational Cost: add. (5032->159), mult. (3914->303), div. (1710->9), fcn. (2844->30), ass. (0->172)
t195 = -2 * pkin(1);
t194 = 2 * pkin(1);
t193 = 2 * pkin(2);
t124 = cos(qJ(3,3));
t192 = t124 * pkin(1);
t126 = cos(qJ(3,2));
t191 = t126 * pkin(1);
t128 = cos(qJ(3,1));
t190 = t128 * pkin(1);
t132 = pkin(2) ^ 2;
t118 = sin(qJ(3,3));
t107 = 0.1e1 / t118;
t133 = 1 / pkin(2);
t134 = 1 / pkin(1);
t148 = t133 * t134;
t130 = xDP(2);
t131 = xDP(1);
t113 = legFrame(3,2);
t97 = -t113 + qJ(2,3);
t91 = qJ(3,3) + t97;
t79 = sin(t91);
t82 = cos(t91);
t52 = t130 * t82 - t131 * t79;
t85 = sin(t97);
t88 = cos(t97);
t37 = -t52 * pkin(2) + (-t130 * t88 + t131 * t85) * pkin(1);
t144 = t37 * t148;
t138 = t107 * t144;
t156 = t107 * t134;
t46 = t52 * t156;
t31 = t46 + t138 / 0.2e1;
t165 = t124 * t31;
t180 = t107 * t52;
t34 = t46 + t138;
t189 = (t34 * t132 + (t165 * t193 + t180) * pkin(1)) * t52;
t120 = sin(qJ(3,2));
t109 = 0.1e1 / t120;
t114 = legFrame(2,2);
t98 = -t114 + qJ(2,2);
t92 = qJ(3,2) + t98;
t80 = sin(t92);
t83 = cos(t92);
t53 = t130 * t83 - t131 * t80;
t86 = sin(t98);
t89 = cos(t98);
t38 = -t53 * pkin(2) + (-t130 * t89 + t131 * t86) * pkin(1);
t143 = t38 * t148;
t137 = t109 * t143;
t154 = t109 * t134;
t47 = t53 * t154;
t32 = t47 + t137 / 0.2e1;
t164 = t126 * t32;
t175 = t109 * t53;
t35 = t47 + t137;
t188 = (t35 * t132 + (t164 * t193 + t175) * pkin(1)) * t53;
t122 = sin(qJ(3,1));
t111 = 0.1e1 / t122;
t115 = legFrame(1,2);
t99 = -t115 + qJ(2,1);
t93 = qJ(3,1) + t99;
t81 = sin(t93);
t84 = cos(t93);
t54 = t130 * t84 - t131 * t81;
t87 = sin(t99);
t90 = cos(t99);
t39 = -t54 * pkin(2) + (-t130 * t90 + t131 * t87) * pkin(1);
t142 = t39 * t148;
t136 = t111 * t142;
t152 = t111 * t134;
t48 = t54 * t152;
t33 = t48 + t136 / 0.2e1;
t163 = t128 * t33;
t170 = t111 * t54;
t36 = t48 + t136;
t187 = (t36 * t132 + (t163 * t193 + t170) * pkin(1)) * t54;
t186 = t34 * t37;
t185 = t35 * t38;
t184 = t36 * t39;
t183 = g(1) * t82 + g(2) * t79;
t182 = g(1) * t83 + g(2) * t80;
t181 = g(1) * t84 + g(2) * t81;
t55 = pkin(1) * t85 + pkin(2) * t79;
t179 = t107 * t55;
t58 = -pkin(1) * t88 - pkin(2) * t82;
t178 = t107 * t58;
t177 = t107 * t79;
t176 = t107 * t82;
t56 = pkin(1) * t86 + pkin(2) * t80;
t174 = t109 * t56;
t59 = -pkin(1) * t89 - pkin(2) * t83;
t173 = t109 * t59;
t172 = t109 * t80;
t171 = t109 * t83;
t57 = pkin(1) * t87 + pkin(2) * t81;
t169 = t111 * t57;
t60 = -pkin(1) * t90 - pkin(2) * t84;
t168 = t111 * t60;
t167 = t111 * t81;
t166 = t111 * t84;
t117 = xDDP(1);
t162 = t55 * t117;
t161 = t56 * t117;
t160 = t57 * t117;
t116 = xDDP(2);
t159 = t58 * t116;
t158 = t59 * t116;
t157 = t60 * t116;
t108 = 0.1e1 / t118 ^ 2;
t135 = 1 / pkin(1) ^ 2;
t155 = t108 * t135;
t110 = 0.1e1 / t120 ^ 2;
t153 = t110 * t135;
t112 = 0.1e1 / t122 ^ 2;
t151 = t112 * t135;
t150 = t116 * t134;
t149 = t117 * t134;
t147 = (pkin(2) + t192) * t186;
t146 = (pkin(2) + t191) * t185;
t145 = (pkin(2) + t190) * t184;
t141 = -g(1) * t79 + g(2) * t82;
t140 = -g(1) * t80 + g(2) * t83;
t139 = -g(1) * t81 + g(2) * t84;
t19 = -t149 * t177 + t150 * t176 + (-(-pkin(2) * t34 - t124 * t180) * t52 + t186) * t155;
t20 = -t149 * t172 + t150 * t171 + (-(-pkin(2) * t35 - t126 * t175) * t53 + t185) * t153;
t21 = -t149 * t167 + t150 * t166 + (-(-t36 * pkin(2) - t128 * t170) * t54 + t184) * t151;
t129 = cos(qJ(2,1));
t127 = cos(qJ(2,2));
t125 = cos(qJ(2,3));
t123 = sin(qJ(2,1));
t121 = sin(qJ(2,2));
t119 = sin(qJ(2,3));
t106 = cos(t115);
t105 = cos(t114);
t104 = cos(t113);
t103 = sin(t115);
t102 = sin(t114);
t101 = sin(t113);
t100 = (xDDP(3) - g(3));
t69 = t106 * g(1) - t103 * g(2);
t68 = t105 * g(1) - t102 * g(2);
t67 = t104 * g(1) - t101 * g(2);
t66 = t103 * g(1) + t106 * g(2);
t65 = t102 * g(1) + t105 * g(2);
t64 = t101 * g(1) + t104 * g(2);
t51 = t54 ^ 2;
t50 = t53 ^ 2;
t49 = t52 ^ 2;
t45 = -t69 * t123 + t66 * t129;
t44 = t66 * t123 + t69 * t129;
t43 = -t68 * t121 + t65 * t127;
t42 = t65 * t121 + t68 * t127;
t41 = -t67 * t119 + t64 * t125;
t40 = t64 * t119 + t67 * t125;
t18 = t51 * t152 + t21 * t190 + t181;
t17 = t49 * t156 + t19 * t192 + t183;
t16 = t50 * t154 + t20 * t191 + t182;
t15 = t128 * t51 * t112 * t134 - t122 * t21 * pkin(1) + t139;
t14 = t126 * t50 * t110 * t134 - t120 * t20 * pkin(1) + t140;
t13 = t124 * t49 * t134 * t108 - t118 * t19 * pkin(1) + t141;
t12 = ((-t145 - t187) * t151 + (t157 + t160) * t152) * t133 + t21;
t11 = ((-t146 - t188) * t153 + (t158 + t161) * t154) * t133 + t20;
t10 = ((-t147 - t189) * t155 + (t159 + t162) * t156) * t133 + t19;
t9 = ((-t187 / 0.2e1 - t145 / 0.2e1) * t151 + (t160 / 0.2e1 + t157 / 0.2e1) * t152) * t133 + t21;
t8 = ((-t188 / 0.2e1 - t146 / 0.2e1) * t153 + (t161 / 0.2e1 + t158 / 0.2e1) * t154) * t133 + t20;
t7 = ((-t189 / 0.2e1 - t147 / 0.2e1) * t155 + (t162 / 0.2e1 + t159 / 0.2e1) * t156) * t133 + t19;
t6 = (t9 * t128 - t33 * t142) * t194 + t181;
t5 = (t8 * t126 - t32 * t143) * t194 + t182;
t4 = (t7 * t124 - t31 * t144) * t194 + t183;
t3 = (t122 * t9 + t136 * t163) * t195 + t139;
t2 = (t120 * t8 + t137 * t164) * t195 + t140;
t1 = (t118 * t7 + t138 * t165) * t195 + t141;
t22 = [0, (-t21 * t167 - t20 * t172 - t19 * t177) * t134, (-t44 * t167 - t42 * t172 - t40 * t177) * t134, (-t45 * t167 - t43 * t172 - t41 * t177) * t134, (-t10 * t177 - t11 * t172 - t12 * t167 + (t10 * t179 + t11 * t174 + t12 * t169) * t133) * t134, (-t4 * t177 - t5 * t172 - t6 * t167 + (t16 * t174 + t18 * t169 + t17 * t179) * t133) * t134, (-t1 * t177 - t2 * t172 - t3 * t167 + (t13 * t179 + t14 * t174 + t15 * t169) * t133) * t134, t117 - g(1); 0, (t21 * t166 + t20 * t171 + t19 * t176) * t134, (t44 * t166 + t42 * t171 + t40 * t176) * t134, (t45 * t166 + t43 * t171 + t41 * t176) * t134, (t10 * t176 + t11 * t171 + t12 * t166 + (t10 * t178 + t11 * t173 + t12 * t168) * t133) * t134, (t4 * t176 + t5 * t171 + t6 * t166 + (t16 * t173 + t18 * t168 + t17 * t178) * t133) * t134, (t1 * t176 + t2 * t171 + t3 * t166 + (t13 * t178 + t14 * t173 + t15 * t168) * t133) * t134, t116 - g(2); 3 * t100, 0, 0, 0, 0, 0, 0, t100;];
tauX_reg  = t22;
