% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR1G3P1A0
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
% koppelP [3x3]
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
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:19
% EndTime: 2020-03-09 21:02:21
% DurationCPUTime: 2.28s
% Computational Cost: add. (1620->267), mult. (3162->538), div. (2328->14), fcn. (3711->18), ass. (0->200)
t109 = cos(qJ(3,1));
t86 = t109 ^ 2;
t213 = (0.2e1 * t86 - 0.1e1) * Ifges(3,4);
t107 = cos(qJ(3,2));
t82 = t107 ^ 2;
t212 = (0.2e1 * t82 - 0.1e1) * Ifges(3,4);
t105 = cos(qJ(3,3));
t78 = t105 ^ 2;
t211 = (0.2e1 * t78 - 0.1e1) * Ifges(3,4);
t210 = 0.2e1 * Ifges(3,4);
t79 = 0.1e1 / t105;
t83 = 0.1e1 / t107;
t87 = 0.1e1 / t109;
t92 = mrSges(2,2) - mrSges(3,3);
t206 = t92 / 0.2e1;
t205 = mrSges(3,1) * g(3);
t100 = sin(qJ(2,3));
t99 = sin(qJ(3,3));
t126 = -mrSges(3,1) * t105 + mrSges(3,2) * t99;
t112 = xDP(3);
t90 = t112 ^ 2;
t168 = 0.1e1 / pkin(2) ^ 2 * t90;
t81 = t79 / t78;
t148 = t81 * t168;
t132 = t99 * t148;
t115 = 0.1e1 / pkin(2);
t160 = t112 * t115;
t142 = t79 * t160;
t106 = cos(qJ(2,3));
t113 = xDP(2);
t114 = xDP(1);
t80 = 0.1e1 / t105 ^ 2;
t158 = t106 * t80 * t99;
t93 = legFrame(3,2);
t68 = sin(t93);
t71 = cos(t93);
t75 = 0.1e1 / t100;
t37 = (-t112 * t158 + (t113 * t68 - t114 * t71) * t79) * t75 * t115;
t177 = t106 * t37;
t201 = mrSges(3,1) * t99;
t169 = t115 * t90;
t34 = t37 ^ 2;
t179 = t105 * t34;
t22 = (-pkin(2) * t179 - t81 * t169) * t75;
t122 = mrSges(2,1) - t126;
t40 = -t100 * t92 + t122 * t106;
t52 = g(1) * t68 + g(2) * t71;
t61 = mrSges(3,2) * t105 + t201;
t67 = mrSges(3,2) * t160;
t172 = t115 * t80;
t178 = t105 * t99;
t7 = ((pkin(2) * t78 * t177 - t100 * t112 * t99) * t37 * t172 + (-t100 * t37 * t178 + t106 * t142) * t81 * t160) * t75;
t74 = m(1) + m(2) + m(3);
t4 = -t40 * t7 - 0.2e1 * (t142 * t201 + t37 * t206 + t67) * t177 + (-t22 - t52) * t74 + (-t61 * t132 - mrSges(2,1) * t34 + t126 * (t80 * t168 + t34)) * t100;
t204 = t75 * t4;
t101 = sin(qJ(3,2));
t102 = sin(qJ(2,2));
t141 = t83 * t160;
t85 = t83 / t82;
t147 = t85 * t168;
t108 = cos(qJ(2,2));
t163 = t101 * t112;
t94 = legFrame(2,2);
t69 = sin(t94);
t72 = cos(t94);
t76 = 0.1e1 / t102;
t84 = 0.1e1 / t107 ^ 2;
t38 = (-t108 * t84 * t163 + (t113 * t69 - t114 * t72) * t83) * t76 * t115;
t175 = t108 * t38;
t189 = mrSges(3,1) * t107;
t191 = mrSges(3,1) * t101;
t35 = t38 ^ 2;
t176 = t107 * t35;
t23 = (-pkin(2) * t176 - t85 * t169) * t76;
t32 = t84 * t168 + t35;
t121 = -mrSges(3,2) * t101 + mrSges(2,1) + t189;
t41 = -t102 * t92 + t121 * t108;
t53 = g(1) * t69 + g(2) * t72;
t62 = mrSges(3,2) * t107 + t191;
t164 = t101 * t107;
t171 = t115 * t84;
t8 = ((pkin(2) * t82 * t175 - t102 * t163) * t38 * t171 + (-t102 * t38 * t164 + t108 * t141) * t85 * t160) * t76;
t203 = t76 * (-t41 * t8 - 0.2e1 * (t141 * t191 + t38 * t206 + t67) * t175 + (-t23 - t53) * t74 + (-t32 * t189 - mrSges(2,1) * t35 + (mrSges(3,2) * t32 - t62 * t147) * t101) * t102);
t103 = sin(qJ(3,1));
t104 = sin(qJ(2,1));
t140 = t87 * t160;
t89 = t87 / t86;
t146 = t89 * t168;
t110 = cos(qJ(2,1));
t161 = t103 * t112;
t95 = legFrame(1,2);
t70 = sin(t95);
t73 = cos(t95);
t77 = 0.1e1 / t104;
t88 = 0.1e1 / t109 ^ 2;
t39 = (-t110 * t88 * t161 + (t113 * t70 - t114 * t73) * t87) * t77 * t115;
t173 = t110 * t39;
t188 = mrSges(3,1) * t109;
t190 = mrSges(3,1) * t103;
t36 = t39 ^ 2;
t174 = t109 * t36;
t24 = (-pkin(2) * t174 - t89 * t169) * t77;
t33 = t88 * t168 + t36;
t120 = -mrSges(3,2) * t103 + mrSges(2,1) + t188;
t42 = -t104 * t92 + t120 * t110;
t54 = g(1) * t70 + g(2) * t73;
t63 = mrSges(3,2) * t109 + t190;
t162 = t103 * t109;
t170 = t115 * t88;
t9 = ((pkin(2) * t86 * t173 - t104 * t161) * t39 * t170 + (-t104 * t39 * t162 + t110 * t140) * t89 * t160) * t77;
t202 = t77 * (-t42 * t9 - 0.2e1 * (t140 * t190 + t39 * t206 + t67) * t173 + (-t24 - t54) * t74 + (-t33 * t188 - mrSges(2,1) * t36 + (mrSges(3,2) * t33 - t63 * t146) * t103) * t104);
t97 = xDDP(2);
t200 = t75 * t97;
t98 = xDDP(1);
t199 = t75 * t98;
t198 = t75 * t99;
t197 = t76 * t97;
t196 = t76 * t98;
t195 = t77 * t97;
t194 = t77 * t98;
t193 = t79 * t99;
t192 = Ifges(3,1) + Ifges(2,3);
t186 = t100 * t61;
t185 = t101 * t76;
t184 = t101 * t83;
t183 = t102 * t62;
t182 = t103 * t77;
t181 = t103 * t87;
t180 = t104 * t63;
t167 = t79 * t115;
t166 = t83 * t115;
t165 = t87 * t115;
t159 = t75 * t193;
t157 = t68 * t167;
t156 = t69 * t166;
t155 = t70 * t165;
t154 = t71 * t167;
t153 = t72 * t166;
t152 = t73 * t165;
t151 = t75 * t167;
t150 = t76 * t166;
t149 = t77 * t165;
t145 = t106 * t172;
t144 = t108 * t171;
t143 = t110 * t170;
t139 = t75 * t158;
t138 = t68 * t151;
t137 = t69 * t150;
t136 = t70 * t149;
t135 = t71 * t151;
t134 = t72 * t150;
t133 = t73 * t149;
t131 = t76 * t144;
t130 = t77 * t143;
t129 = Ifges(3,5) * t160 / 0.2e1;
t128 = -Ifges(3,6) * t160 / 0.2e1;
t127 = t115 * t139;
t55 = g(1) * t71 - g(2) * t68;
t125 = t100 * t52 + t106 * t55;
t56 = g(1) * t72 - g(2) * t69;
t124 = t102 * t53 + t108 * t56;
t57 = g(1) * t73 - g(2) * t70;
t123 = t104 * t54 + t110 * t57;
t111 = mrSges(3,2) * g(3);
t96 = xDDP(3);
t91 = Ifges(3,1) - Ifges(3,2);
t60 = Ifges(3,5) * t103 + Ifges(3,6) * t109;
t59 = Ifges(3,5) * t101 + Ifges(3,6) * t107;
t58 = Ifges(3,5) * t99 + Ifges(3,6) * t105;
t51 = t104 * t70 + t110 * t73;
t50 = t104 * t73 - t110 * t70;
t49 = t102 * t69 + t108 * t72;
t48 = t102 * t72 - t108 * t69;
t47 = t100 * t68 + t106 * t71;
t46 = t100 * t71 - t106 * t68;
t45 = t162 * t210 - t86 * t91 + t192;
t44 = t164 * t210 - t82 * t91 + t192;
t43 = t178 * t210 - t78 * t91 + t192;
t30 = (-t42 * t152 + t51 * t74) * t77;
t29 = (t42 * t155 + t50 * t74) * t77;
t28 = (-t41 * t153 + t49 * t74) * t76;
t27 = (t41 * t156 + t48 * t74) * t76;
t26 = (-t40 * t154 + t47 * t74) * t75;
t25 = (t40 * t157 + t46 * t74) * t75;
t21 = -t165 * t180 + (-t42 * t143 + t74 * t87) * t182;
t20 = -t166 * t183 + (-t41 * t144 + t74 * t83) * t185;
t19 = t74 * t159 + (-t40 * t139 - t79 * t186) * t115;
t18 = (-t45 * t152 + t42 * t51) * t77;
t17 = (t45 * t155 + t42 * t50) * t77;
t16 = (-t44 * t153 + t41 * t49) * t76;
t15 = (t44 * t156 + t41 * t48) * t76;
t14 = (-t43 * t154 + t40 * t47) * t75;
t13 = (t43 * t157 + t40 * t46) * t75;
t12 = t60 * t165 + (-t45 * t143 + t42 * t87) * t182;
t11 = t59 * t166 + (-t44 * t144 + t41 * t83) * t185;
t10 = t40 * t159 + (-t43 * t139 + t58 * t79) * t115;
t3 = -t42 * t24 - t45 * t9 + t60 * t103 * t146 + 0.2e1 * ((t39 * t91 * t103 + t87 * t129) * t109 + t128 * t181 + t39 * t213) * t140 + (-t120 * t54 + t57 * t92) * t110 + (t120 * t57 + t54 * t92) * t104;
t2 = -t41 * t23 - t44 * t8 + t59 * t101 * t147 + 0.2e1 * ((t38 * t91 * t101 + t83 * t129) * t107 + t128 * t184 + t38 * t212) * t141 + (-t121 * t53 + t56 * t92) * t108 + (t121 * t56 + t53 * t92) * t102;
t1 = -t40 * t22 - t43 * t7 + t58 * t132 + 0.2e1 * ((t37 * t91 * t99 + t79 * t129) * t105 + t128 * t193 + t37 * t211) * t142 + (-t122 * t52 + t55 * t92) * t106 + (t122 * t55 + t52 * t92) * t100;
t5 = [(-t18 * t152 + t30 * t51) * t194 + (t18 * t155 + t30 * t50) * t195 + t51 * t202 - t3 * t133 + (-t16 * t153 + t28 * t49) * t196 + (t16 * t156 + t28 * t48) * t197 + t49 * t203 - t2 * t134 + (-t14 * t154 + t26 * t47) * t199 + (t14 * t157 + t26 * t46) * t200 + t47 * t204 - t1 * t135 + (-g(1) + t98) * m(4) + ((-t60 * t133 - t51 * t63) * t165 + (-t18 * t143 + t30 * t87) * t182 + (-t59 * t134 - t49 * t62) * t166 + (-t16 * t144 + t28 * t83) * t185 + (-t58 * t135 - t47 * t61) * t167 + (-t14 * t145 + t26 * t79) * t198) * t96; (-t17 * t152 + t29 * t51) * t194 + (t17 * t155 + t29 * t50) * t195 + t50 * t202 + t3 * t136 + (-t15 * t153 + t27 * t49) * t196 + (t15 * t156 + t27 * t48) * t197 + t48 * t203 + t2 * t137 + (-t13 * t154 + t25 * t47) * t199 + (t13 * t157 + t25 * t46) * t200 + t46 * t204 + t1 * t138 + (-g(2) + t97) * m(4) + ((t60 * t136 - t50 * t63) * t165 + (-t17 * t143 + t29 * t87) * t182 + (t59 * t137 - t48 * t62) * t166 + (-t15 * t144 + t27 * t83) * t185 + (t58 * t138 - t46 * t61) * t167 + (-t13 * t145 + t25 * t79) * t198) * t96; (-t12 * t152 + t21 * t51) * t194 + (t12 * t155 + t21 * t50) * t195 + t181 * t202 - t103 * t3 * t130 + (t24 * t180 - t60 * t9 + (t123 * mrSges(3,2) - t205) * t109 + (t123 * mrSges(3,1) + Ifges(3,3) * t146 - t91 * t174 + t111) * t103 - t36 * t213) * t165 + (-t11 * t153 + t20 * t49) * t196 + (t11 * t156 + t20 * t48) * t197 + t184 * t203 - t101 * t2 * t131 + (t23 * t183 - t59 * t8 + (t124 * mrSges(3,2) - t205) * t107 + (t124 * mrSges(3,1) + Ifges(3,3) * t147 - t91 * t176 + t111) * t101 - t35 * t212) * t166 + (-t10 * t154 + t19 * t47) * t199 + (t10 * t157 + t19 * t46) * t200 + t4 * t159 - t1 * t127 + (t22 * t186 - t58 * t7 + (t125 * mrSges(3,2) - t205) * t105 + (t125 * mrSges(3,1) + Ifges(3,3) * t148 - t91 * t179 + t111) * t99 - t34 * t211) * t167 - g(3) * m(4) + ((t21 * t87 * t77 - t12 * t130 + (-t60 * t130 - t63 * t87) * t165) * t103 + (t20 * t83 * t76 - t11 * t131 + (-t59 * t131 - t62 * t83) * t166) * t101 + t19 * t159 + (-t10 * t139 + (-t58 * t127 + (Ifges(3,3) * t115 - t61 * t99) * t79) * t79) * t115 + m(4) + (t83 ^ 2 + t87 ^ 2) * Ifges(3,3) * t115 ^ 2) * t96;];
tauX  = t5;
