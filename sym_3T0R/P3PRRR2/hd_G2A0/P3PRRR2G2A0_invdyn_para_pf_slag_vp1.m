% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR2G2A0
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
%   pkin=[a3,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR2G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:29
% EndTime: 2020-03-09 21:21:31
% DurationCPUTime: 1.80s
% Computational Cost: add. (3645->233), mult. (6804->415), div. (2988->5), fcn. (6144->24), ass. (0->193)
t132 = xDP(2);
t111 = legFrame(1,2);
t98 = cos(t111);
t182 = t98 * t132;
t133 = xDP(1);
t95 = sin(t111);
t185 = t95 * t133;
t216 = t182 + t185;
t110 = legFrame(2,2);
t97 = cos(t110);
t183 = t97 * t132;
t94 = sin(t110);
t186 = t94 * t133;
t215 = t183 + t186;
t109 = legFrame(3,2);
t96 = cos(t109);
t184 = t96 * t132;
t93 = sin(t109);
t187 = t93 * t133;
t214 = t184 + t187;
t213 = m(3) * pkin(1);
t121 = cos(qJ(3,3));
t212 = pkin(1) * t121;
t123 = cos(qJ(3,2));
t211 = pkin(1) * t123;
t125 = cos(qJ(3,1));
t210 = pkin(1) * t125;
t137 = 0.1e1 / pkin(2);
t122 = cos(qJ(2,3));
t115 = sin(qJ(3,3));
t116 = sin(qJ(2,3));
t178 = t115 * t116;
t64 = (pkin(2) * t121 + pkin(1)) * t122 - pkin(2) * t178;
t197 = t137 * t64;
t148 = t214 * t197;
t103 = 0.1e1 / t115;
t139 = 0.1e1 / pkin(1);
t181 = t103 * t139;
t131 = xDP(3);
t106 = qJ(2,3) + qJ(3,3);
t87 = sin(t106);
t190 = t87 * t131;
t67 = t121 * t122 - t178;
t206 = t214 * t67 * t181;
t175 = t131 * t139;
t166 = t137 * t175;
t76 = pkin(1) * t116 + pkin(2) * t87;
t55 = t76 * t103 * t166;
t16 = t55 + (-t148 - t190) * t181 + t206;
t25 = -t148 * t181 + t55;
t209 = t16 * t25;
t124 = cos(qJ(2,2));
t117 = sin(qJ(3,2));
t118 = sin(qJ(2,2));
t177 = t117 * t118;
t65 = (pkin(2) * t123 + pkin(1)) * t124 - pkin(2) * t177;
t196 = t137 * t65;
t147 = t215 * t196;
t104 = 0.1e1 / t117;
t180 = t104 * t139;
t107 = qJ(2,2) + qJ(3,2);
t88 = sin(t107);
t189 = t88 * t131;
t68 = t123 * t124 - t177;
t205 = t215 * t68 * t180;
t77 = pkin(1) * t118 + pkin(2) * t88;
t56 = t77 * t104 * t166;
t17 = t56 + (-t147 - t189) * t180 + t205;
t26 = -t147 * t180 + t56;
t208 = t17 * t26;
t126 = cos(qJ(2,1));
t119 = sin(qJ(3,1));
t120 = sin(qJ(2,1));
t176 = t119 * t120;
t66 = (pkin(2) * t125 + pkin(1)) * t126 - pkin(2) * t176;
t195 = t137 * t66;
t146 = t216 * t195;
t105 = 0.1e1 / t119;
t179 = t105 * t139;
t108 = qJ(2,1) + qJ(3,1);
t89 = sin(t108);
t188 = t89 * t131;
t69 = t125 * t126 - t176;
t204 = t216 * t69 * t179;
t78 = pkin(1) * t120 + pkin(2) * t89;
t57 = t78 * t105 * t166;
t18 = t57 + (-t146 - t188) * t179 + t204;
t27 = -t146 * t179 + t57;
t207 = t18 * t27;
t113 = xDDP(2);
t203 = t113 * t96;
t202 = t113 * t97;
t201 = t113 * t98;
t114 = xDDP(1);
t200 = t114 * t93;
t199 = t114 * t94;
t198 = t114 * t95;
t194 = t137 * t76;
t193 = t137 * t77;
t192 = t137 * t78;
t99 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t82 = m(3) * t99 + Icges(3,3);
t191 = t137 * t82;
t174 = 0.2e1 * pkin(1) * pkin(2);
t173 = rSges(3,2) * t213;
t172 = rSges(3,1) * t212;
t171 = rSges(3,1) * t211;
t170 = rSges(3,1) * t210;
t165 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3);
t127 = rSges(3,2) * g(3);
t129 = rSges(3,1) * g(3);
t70 = g(1) * t93 + g(2) * t96;
t164 = (rSges(3,1) * t70 - t127) * t87 + (rSges(3,2) * t70 + t129) * cos(t106);
t71 = g(1) * t94 + g(2) * t97;
t163 = (rSges(3,1) * t71 - t127) * t88 + (rSges(3,2) * t71 + t129) * cos(t107);
t72 = g(1) * t95 + g(2) * t98;
t162 = (rSges(3,1) * t72 - t127) * t89 + (rSges(3,2) * t72 + t129) * cos(t108);
t37 = -t103 * t175 * t87 + t206;
t10 = ((-pkin(2) * t16 - t37 * t212) * t37 - pkin(2) * t209) * t181;
t102 = g(3) * t213;
t128 = rSges(2,2) * g(3);
t13 = t55 / 0.2e1 + (-t190 + (-t187 / 0.2e1 - t184 / 0.2e1) * t197) * t181 + t206;
t130 = rSges(2,1) * g(3);
t83 = t115 * t173;
t138 = pkin(1) ^ 2;
t86 = t138 + t99;
t43 = -0.2e1 * t83 + (t86 + 0.2e1 * t172) * m(3) + t165;
t52 = Icges(3,3) - t83 + (t99 + t172) * m(3);
t136 = pkin(2) ^ 2;
t7 = ((t121 * t13 * t174 + t136 * t16 + t138 * t37) * t137 * t37 + (pkin(2) + t212) * t209) * t181;
t79 = rSges(3,1) * t115 + rSges(3,2) * t121;
t1 = -t43 * t10 - t52 * t7 + (t102 + m(2) * (rSges(2,2) * t70 + t130)) * t122 - m(2) * (-rSges(2,1) * t70 + t128) * t116 + ((-0.2e1 * t13 * t25 * t79 + t116 * t70) * pkin(1) + t164) * m(3);
t4 = -t52 * t10 - t82 * t7 + (pkin(1) * t37 ^ 2 * t79 + t164) * m(3);
t161 = t67 * t1 - t4 * t197;
t38 = -t104 * t175 * t88 + t205;
t11 = ((-pkin(2) * t17 - t38 * t211) * t38 - pkin(2) * t208) * t180;
t14 = t56 / 0.2e1 + (-t189 + (-t186 / 0.2e1 - t183 / 0.2e1) * t196) * t180 + t205;
t84 = t117 * t173;
t44 = -0.2e1 * t84 + (t86 + 0.2e1 * t171) * m(3) + t165;
t53 = Icges(3,3) - t84 + (t99 + t171) * m(3);
t8 = ((t123 * t14 * t174 + t136 * t17 + t138 * t38) * t137 * t38 + (pkin(2) + t211) * t208) * t180;
t80 = rSges(3,1) * t117 + rSges(3,2) * t123;
t2 = -t44 * t11 - t53 * t8 + (t102 + m(2) * (rSges(2,2) * t71 + t130)) * t124 - m(2) * (-rSges(2,1) * t71 + t128) * t118 + ((-0.2e1 * t14 * t26 * t80 + t118 * t71) * pkin(1) + t163) * m(3);
t5 = -t53 * t11 - t82 * t8 + (pkin(1) * t38 ^ 2 * t80 + t163) * m(3);
t160 = -t5 * t196 + t68 * t2;
t39 = -t105 * t175 * t89 + t204;
t12 = ((-pkin(2) * t18 - t39 * t210) * t39 - pkin(2) * t207) * t179;
t15 = t57 / 0.2e1 + (-t188 + (-t185 / 0.2e1 - t182 / 0.2e1) * t195) * t179 + t204;
t85 = t119 * t173;
t45 = -0.2e1 * t85 + (t86 + 0.2e1 * t170) * m(3) + t165;
t54 = Icges(3,3) - t85 + (t99 + t170) * m(3);
t81 = rSges(3,1) * t119 + rSges(3,2) * t125;
t9 = ((t125 * t15 * t174 + t136 * t18 + t138 * t39) * t137 * t39 + (pkin(2) + t210) * t207) * t179;
t3 = -t45 * t12 - t54 * t9 + (t102 + m(2) * (rSges(2,2) * t72 + t130)) * t126 - m(2) * (-rSges(2,1) * t72 + t128) * t120 + ((-0.2e1 * t15 * t27 * t81 + t120 * t72) * pkin(1) + t162) * m(3);
t6 = -t54 * t12 - t82 * t9 + (pkin(1) * t39 ^ 2 * t81 + t162) * m(3);
t159 = -t6 * t195 + t69 * t3;
t145 = (-t52 * t197 + t43 * t67) * t181;
t19 = t93 * t145;
t144 = (-t64 * t191 + t52 * t67) * t181;
t28 = t93 * t144;
t158 = t19 * t67 - t28 * t197;
t143 = (-t53 * t196 + t44 * t68) * t180;
t20 = t94 * t143;
t142 = (-t65 * t191 + t53 * t68) * t180;
t29 = t94 * t142;
t157 = -t29 * t196 + t20 * t68;
t141 = (-t54 * t195 + t45 * t69) * t179;
t21 = t95 * t141;
t140 = (-t66 * t191 + t54 * t69) * t179;
t30 = t95 * t140;
t156 = -t30 * t195 + t21 * t69;
t22 = t96 * t145;
t31 = t96 * t144;
t155 = -t31 * t197 + t22 * t67;
t23 = t97 * t143;
t32 = t97 * t142;
t154 = -t32 * t196 + t23 * t68;
t24 = t98 * t141;
t33 = t98 * t140;
t153 = -t33 * t195 + t24 * t69;
t149 = -t93 * t96 - t94 * t97 - t95 * t98;
t112 = xDDP(3);
t101 = m(1) + m(2) + m(3);
t75 = g(1) * t98 - g(2) * t95;
t74 = g(1) * t97 - g(2) * t94;
t73 = g(1) * t96 - g(2) * t93;
t42 = (t78 * t191 - t54 * t89) * t179;
t41 = (t77 * t191 - t53 * t88) * t180;
t40 = (t76 * t191 - t52 * t87) * t181;
t36 = (t54 * t192 - t45 * t89) * t179;
t35 = (t53 * t193 - t44 * t88) * t180;
t34 = (t52 * t194 - t43 * t87) * t181;
t46 = [(-g(1) + t114) * m(4) + (-t96 * t73 - t97 * t74 - t98 * t75 + (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) * t114 + t149 * t113) * t101 + (((t30 * t192 - t21 * t89) * t112 + t156 * t201 + (t114 * t156 + t159) * t95) * t105 + ((t29 * t193 - t20 * t88) * t112 + t157 * t202 + (t114 * t157 + t160) * t94) * t104 + ((-t19 * t87 + t28 * t194) * t112 + t158 * t203 + (t114 * t158 + t161) * t93) * t103) * t139; (-g(2) + t113) * m(4) + (t93 * t73 + t94 * t74 + t95 * t75 + t149 * t114 + (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) * t113) * t101 + (((t33 * t192 - t24 * t89) * t112 + t153 * t198 + (t113 * t153 + t159) * t98) * t105 + ((t32 * t193 - t23 * t88) * t112 + t154 * t199 + (t113 * t154 + t160) * t97) * t104 + ((t31 * t194 - t22 * t87) * t112 + t155 * t200 + (t113 * t155 + t161) * t96) * t103) * t139; (-g(3) + t112) * m(4) + (((t42 * t192 - t36 * t89) * t112 - t89 * t3 + t6 * t192 + (t201 + t198) * (-t42 * t195 + t36 * t69)) * t105 + ((t41 * t193 - t35 * t88) * t112 - t88 * t2 + t5 * t193 + (t202 + t199) * (-t41 * t196 + t35 * t68)) * t104 + ((t40 * t194 - t34 * t87) * t112 - t87 * t1 + t4 * t194 + (t203 + t200) * (-t40 * t197 + t34 * t67)) * t103) * t139;];
tauX  = t46;
