% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR1G1A0
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:04
% EndTime: 2020-03-09 20:34:06
% DurationCPUTime: 1.59s
% Computational Cost: add. (1782->248), mult. (4060->493), div. (2172->10), fcn. (4935->18), ass. (0->221)
t132 = cos(qJ(3,1));
t111 = t132 ^ 2;
t242 = (0.2e1 * t111 - 0.1e1) * Ifges(3,4);
t130 = cos(qJ(3,2));
t108 = t130 ^ 2;
t241 = (0.2e1 * t108 - 0.1e1) * Ifges(3,4);
t128 = cos(qJ(3,3));
t105 = t128 ^ 2;
t240 = (0.2e1 * t105 - 0.1e1) * Ifges(3,4);
t239 = 0.2e1 * Ifges(3,4);
t235 = Ifges(3,5) / 0.2e1;
t234 = -Ifges(3,6) / 0.2e1;
t118 = mrSges(2,2) - mrSges(3,3);
t233 = t118 / 0.2e1;
t101 = m(1) + m(2) + m(3);
t232 = g(3) * t101;
t231 = Ifges(3,1) + Ifges(2,3);
t136 = 0.1e1 / pkin(2);
t229 = Ifges(3,3) * t136;
t106 = 0.1e1 / t128;
t134 = xDP(2);
t135 = xDP(1);
t198 = t106 * t136;
t114 = legFrame(3,3);
t94 = sin(t114);
t97 = cos(t114);
t58 = (-t134 * t97 + t135 * t94) * t198;
t55 = t58 ^ 2;
t228 = t106 * t55;
t123 = sin(qJ(2,3));
t129 = cos(qJ(2,3));
t122 = sin(qJ(3,3));
t145 = -mrSges(3,1) * t128 + mrSges(3,2) * t122;
t142 = mrSges(2,1) - t145;
t61 = -t123 * t118 + t142 * t129;
t227 = t106 * t61;
t226 = t106 * t94;
t225 = t106 * t97;
t109 = 0.1e1 / t130;
t196 = t109 * t136;
t115 = legFrame(2,3);
t95 = sin(t115);
t98 = cos(t115);
t59 = (-t134 * t98 + t135 * t95) * t196;
t56 = t59 ^ 2;
t224 = t109 * t56;
t125 = sin(qJ(2,2));
t131 = cos(qJ(2,2));
t124 = sin(qJ(3,2));
t144 = -mrSges(3,1) * t130 + mrSges(3,2) * t124;
t141 = mrSges(2,1) - t144;
t62 = -t125 * t118 + t141 * t131;
t223 = t109 * t62;
t222 = t109 * t95;
t221 = t109 * t98;
t112 = 0.1e1 / t132;
t194 = t112 * t136;
t116 = legFrame(1,3);
t96 = sin(t116);
t99 = cos(t116);
t60 = (-t134 * t99 + t135 * t96) * t194;
t57 = t60 ^ 2;
t220 = t112 * t57;
t127 = sin(qJ(2,1));
t133 = cos(qJ(2,1));
t126 = sin(qJ(3,1));
t143 = -mrSges(3,1) * t132 + mrSges(3,2) * t126;
t140 = mrSges(2,1) - t143;
t63 = -t127 * t118 + t140 * t133;
t219 = t112 * t63;
t218 = t112 * t96;
t217 = t112 * t99;
t216 = t122 * t58;
t88 = mrSges(3,1) * t122 + mrSges(3,2) * t128;
t215 = t123 * t88;
t214 = t124 * t59;
t89 = mrSges(3,1) * t124 + mrSges(3,2) * t130;
t213 = t125 * t89;
t212 = t126 * t60;
t90 = mrSges(3,1) * t126 + mrSges(3,2) * t132;
t211 = t127 * t90;
t102 = 0.1e1 / t123;
t107 = 0.1e1 / t128 ^ 2;
t197 = t107 * t136;
t157 = t102 * t197;
t191 = t122 * t129;
t64 = -t128 * t97 - t94 * t191;
t68 = -t128 * t94 + t97 * t191;
t52 = (t134 * t68 + t135 * t64) * t157;
t49 = t52 ^ 2;
t210 = t128 * t49;
t103 = 0.1e1 / t125;
t110 = 0.1e1 / t130 ^ 2;
t195 = t110 * t136;
t156 = t103 * t195;
t189 = t124 * t131;
t65 = -t130 * t98 - t95 * t189;
t70 = -t130 * t95 + t98 * t189;
t53 = (t134 * t70 + t135 * t65) * t156;
t50 = t53 ^ 2;
t209 = t130 * t50;
t104 = 0.1e1 / t127;
t113 = 0.1e1 / t132 ^ 2;
t193 = t113 * t136;
t155 = t104 * t193;
t187 = t126 * t133;
t66 = -t132 * t99 - t96 * t187;
t72 = -t132 * t96 + t99 * t187;
t54 = (t134 * t72 + t135 * t66) * t155;
t51 = t54 ^ 2;
t208 = t132 * t51;
t207 = t101 * t106;
t206 = t101 * t109;
t205 = t101 * t112;
t204 = t102 * t106;
t203 = t102 * t107;
t202 = t103 * t109;
t201 = t103 * t110;
t200 = t104 * t112;
t199 = t104 * t113;
t192 = t122 * t128;
t190 = t124 * t130;
t188 = t126 * t132;
t186 = t128 * t129;
t185 = t130 * t131;
t184 = t132 * t133;
t169 = t106 * t215;
t151 = t136 * t169;
t167 = t61 * t197;
t67 = t122 * t94 + t97 * t186;
t31 = -t94 * t151 + (t64 * t167 + t67 * t207) * t102;
t165 = t109 * t213;
t150 = t136 * t165;
t163 = t62 * t195;
t69 = t124 * t95 + t98 * t185;
t32 = -t95 * t150 + (t65 * t163 + t69 * t206) * t103;
t161 = t112 * t211;
t149 = t136 * t161;
t159 = t63 * t193;
t71 = t126 * t96 + t99 * t184;
t33 = -t96 * t149 + (t66 * t159 + t71 * t205) * t104;
t183 = t33 + t32 + t31;
t73 = -t122 * t97 + t94 * t186;
t34 = t97 * t151 + (t68 * t167 + t73 * t207) * t102;
t74 = -t124 * t98 + t95 * t185;
t35 = t98 * t150 + (t70 * t163 + t74 * t206) * t103;
t75 = -t126 * t99 + t96 * t184;
t36 = t99 * t149 + (t72 * t159 + t75 * t205) * t104;
t182 = t36 + t35 + t34;
t181 = t67 * t204;
t180 = t73 * t204;
t179 = t64 * t203;
t178 = t68 * t203;
t177 = t69 * t202;
t176 = t74 * t202;
t175 = t65 * t201;
t174 = t70 * t201;
t173 = t71 * t200;
t172 = t75 * t200;
t171 = t66 * t199;
t170 = t72 * t199;
t85 = Ifges(3,5) * t122 + Ifges(3,6) * t128;
t168 = t85 * t198;
t117 = Ifges(3,1) - Ifges(3,2);
t76 = -t105 * t117 + t192 * t239 + t231;
t166 = t76 * t197;
t86 = Ifges(3,5) * t124 + Ifges(3,6) * t130;
t164 = t86 * t196;
t77 = -t108 * t117 + t190 * t239 + t231;
t162 = t77 * t195;
t87 = Ifges(3,5) * t126 + Ifges(3,6) * t132;
t160 = t87 * t194;
t78 = -t111 * t117 + t188 * t239 + t231;
t158 = t78 * t193;
t154 = t85 * t157;
t153 = t86 * t156;
t152 = t87 * t155;
t82 = g(1) * t97 + g(2) * t94;
t148 = g(3) * t123 + t129 * t82;
t83 = g(1) * t98 + g(2) * t95;
t147 = g(3) * t125 + t131 * t83;
t84 = g(1) * t99 + g(2) * t96;
t146 = g(3) * t127 + t133 * t84;
t15 = ((-t123 * t216 + t52 * t186) * t106 * t52 + (-t123 * t52 * t192 + t129 * t58) * t107 * t58) * t102;
t37 = (-t210 - t228) * t102 * pkin(2);
t139 = -t122 * t55 * t169 - t61 * t15 + (-mrSges(2,1) * t49 + t145 * (t49 + t55)) * t123 - 0.2e1 * t52 * t129 * (t52 * t233 + t88 * t58) - t101 * t37;
t13 = ((-t125 * t214 + t53 * t185) * t109 * t53 + (-t125 * t53 * t190 + t131 * t59) * t110 * t59) * t103;
t38 = (-t209 - t224) * t103 * pkin(2);
t138 = -t124 * t56 * t165 - t62 * t13 + (-mrSges(2,1) * t50 + t144 * (t50 + t56)) * t125 - 0.2e1 * t53 * t131 * (t53 * t233 + t89 * t59) - t101 * t38;
t14 = ((-t127 * t212 + t54 * t184) * t112 * t54 + (-t127 * t54 * t188 + t133 * t60) * t113 * t60) * t104;
t39 = (-t208 - t220) * t104 * pkin(2);
t137 = -t126 * t57 * t161 - t63 * t14 + (-mrSges(2,1) * t51 + t143 * (t51 + t57)) * t127 - 0.2e1 * t54 * t133 * (t54 * t233 + t90 * t60) - t101 * t39;
t121 = xDDP(1);
t120 = xDDP(2);
t119 = xDDP(3);
t100 = g(3) * t118;
t81 = -g(1) * t96 + g(2) * t99;
t80 = -g(1) * t95 + g(2) * t98;
t79 = -g(1) * t94 + g(2) * t97;
t48 = t72 * t152 + (-t99 * t229 - t75 * t90) * t112;
t47 = t70 * t153 + (-t98 * t229 - t74 * t89) * t109;
t46 = t68 * t154 + (-t97 * t229 - t73 * t88) * t106;
t45 = t66 * t152 + (t96 * t229 - t71 * t90) * t112;
t44 = t65 * t153 + (t95 * t229 - t69 * t89) * t109;
t43 = t64 * t154 + (t94 * t229 - t67 * t88) * t106;
t27 = -t99 * t160 + (t72 * t158 + t75 * t219) * t104;
t26 = -t98 * t164 + (t70 * t162 + t74 * t223) * t103;
t25 = -t97 * t168 + (t68 * t166 + t73 * t227) * t102;
t24 = t96 * t160 + (t66 * t158 + t71 * t219) * t104;
t23 = t95 * t164 + (t65 * t162 + t69 * t223) * t103;
t22 = t94 * t168 + (t64 * t166 + t67 * t227) * t102;
t9 = t39 * t211 - t87 * t14 + (mrSges(3,1) * t81 + t146 * mrSges(3,2)) * t132 + (t146 * mrSges(3,1) - mrSges(3,2) * t81 + Ifges(3,3) * t220 - t117 * t208) * t126 - t51 * t242;
t8 = t38 * t213 - t86 * t13 + (mrSges(3,1) * t80 + t147 * mrSges(3,2)) * t130 + (t147 * mrSges(3,1) - mrSges(3,2) * t80 + Ifges(3,3) * t224 - t117 * t209) * t124 - t50 * t241;
t7 = t37 * t215 - t85 * t15 + (mrSges(3,1) * t79 + t148 * mrSges(3,2)) * t128 + (t148 * mrSges(3,1) - mrSges(3,2) * t79 + Ifges(3,3) * t228 - t117 * t210) * t122 - t49 * t240;
t6 = -t63 * t39 - t78 * t14 + t87 * t126 * t220 + 0.2e1 * ((t54 * t117 * t126 + t60 * t235) * t132 + t212 * t234 + t54 * t242) * t60 + (-t140 * g(3) + t118 * t84) * t133 + t127 * (t140 * t84 + t100);
t5 = -t62 * t38 - t77 * t13 + t86 * t124 * t224 + 0.2e1 * ((t53 * t117 * t124 + t59 * t235) * t130 + t214 * t234 + t53 * t241) * t59 + (-t141 * g(3) + t118 * t83) * t131 + t125 * (t141 * t83 + t100);
t4 = -t61 * t37 - t76 * t15 + t85 * t122 * t228 + 0.2e1 * ((t52 * t117 * t122 + t58 * t235) * t128 + t216 * t234 + t52 * t240) * t58 + (-t142 * g(3) + t118 * t82) * t129 + t123 * (t142 * t82 + t100);
t3 = t137 - t232;
t2 = t138 - t232;
t1 = t139 - t232;
t10 = [t1 * t181 + t2 * t177 + t3 * t173 - g(1) * m(4) + (t33 * t172 + t32 * t176 + t31 * t180) * t120 + t183 * t119 + (t33 * t173 + t32 * t177 + t31 * t181 + m(4)) * t121 + (t4 * t179 + t5 * t175 + t6 * t171 + t7 * t226 + t8 * t222 + t9 * t218 + (t24 * t171 + t23 * t175 + t22 * t179 + t45 * t218 + t44 * t222 + t43 * t226) * t121 + (t24 * t170 + t23 * t174 + t22 * t178 - t45 * t217 - t44 * t221 - t43 * t225) * t120) * t136; t1 * t180 + t2 * t176 + t3 * t172 - g(2) * m(4) + (t36 * t173 + t35 * t177 + t34 * t181) * t121 + t182 * t119 + (t36 * t172 + t35 * t176 + t34 * t180 + m(4)) * t120 + (t4 * t178 + t5 * t174 + t6 * t170 - t7 * t225 - t8 * t221 - t9 * t217 + (t27 * t171 + t26 * t175 + t25 * t179 + t48 * t218 + t47 * t222 + t46 * t226) * t121 + (t27 * t170 + t26 * t174 + t25 * t178 - t48 * t217 - t47 * t221 - t46 * t225) * t120) * t136; t182 * t120 + t183 * t121 + t137 + t138 + t139 + (t119 - g(3)) * (0.3e1 * t101 + m(4));];
tauX  = t10;
