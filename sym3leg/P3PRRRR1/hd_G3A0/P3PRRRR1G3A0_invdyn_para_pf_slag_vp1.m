% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR1G3A0
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
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:14
% EndTime: 2020-03-09 21:02:17
% DurationCPUTime: 2.30s
% Computational Cost: add. (2100->298), mult. (4806->602), div. (2328->14), fcn. (3711->24), ass. (0->213)
t78 = (m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4));
t230 = 2 * t78;
t114 = cos(qJ(3,3));
t90 = 0.1e1 / t114;
t116 = cos(qJ(3,2));
t94 = 0.1e1 / t116;
t118 = cos(qJ(3,1));
t98 = 0.1e1 / t118;
t130 = rSges(3,2) ^ 2;
t131 = rSges(3,1) ^ 2;
t68 = (-t130 + t131) * m(3) - Icges(3,1) + Icges(3,2);
t229 = t68 / 0.2e1;
t223 = rSges(3,3) * m(3);
t76 = rSges(3,2) * t223 - Icges(3,6);
t228 = -t76 / 0.4e1;
t77 = rSges(3,1) * t223 - Icges(3,5);
t227 = t77 / 0.4e1;
t123 = m(2) * rSges(2,1);
t108 = sin(qJ(3,3));
t216 = rSges(3,2) * t114;
t61 = rSges(3,1) * t108 + t216;
t226 = m(3) * t61;
t110 = sin(qJ(3,2));
t215 = rSges(3,2) * t116;
t62 = rSges(3,1) * t110 + t215;
t225 = m(3) * t62;
t112 = sin(qJ(3,1));
t214 = rSges(3,2) * t118;
t63 = rSges(3,1) * t112 + t214;
t224 = m(3) * t63;
t109 = sin(qJ(2,3));
t124 = xDP(3);
t132 = 0.1e1 / pkin(2);
t183 = t124 * t132;
t157 = 0.2e1 * m(3) * t183;
t144 = rSges(3,1) * t157;
t101 = t124 ^ 2;
t190 = t101 / pkin(2) ^ 2;
t89 = t114 ^ 2;
t92 = t90 / t89;
t169 = t92 * t190;
t194 = t90 * t108;
t115 = cos(qJ(2,3));
t125 = xDP(2);
t126 = xDP(1);
t188 = t108 * t124;
t102 = legFrame(3,2);
t79 = sin(t102);
t82 = cos(t102);
t86 = 0.1e1 / t109;
t91 = 0.1e1 / t114 ^ 2;
t40 = (-t115 * t91 * t188 + (t125 * t79 - t126 * t82) * t90) * t86 * t132;
t204 = t115 * t40;
t219 = rSges(3,1) * t114;
t191 = t101 * t132;
t37 = t40 ^ 2;
t28 = (-pkin(2) * t114 * t37 - t92 * t191) * t86;
t31 = t91 * t190 + t37;
t143 = -rSges(3,2) * t108 + t219;
t75 = m(2) * rSges(2,2) - t223;
t43 = (t143 * m(3) + t123) * t115 - t109 * t75;
t55 = g(1) * t79 + g(2) * t82;
t74 = rSges(3,2) * t157;
t85 = m(1) + m(2) + m(3);
t161 = t90 * t183;
t189 = t108 * t114;
t200 = t132 * t91;
t9 = ((pkin(2) * t89 * t204 - t109 * t188) * t40 * t200 + (-t109 * t40 * t189 + t115 * t161) * t92 * t183) * t86;
t222 = (-t43 * t9 - (t144 * t194 + t40 * t75 + t74) * t204 + (-t28 - t55) * t85 + (-t37 * t123 + (-t31 * t219 + (rSges(3,2) * t31 - t61 * t169) * t108) * m(3)) * t109) * t86;
t111 = sin(qJ(2,2));
t93 = t116 ^ 2;
t96 = t94 / t93;
t168 = t96 * t190;
t193 = t94 * t110;
t117 = cos(qJ(2,2));
t186 = t110 * t124;
t103 = legFrame(2,2);
t80 = sin(t103);
t83 = cos(t103);
t87 = 0.1e1 / t111;
t95 = 0.1e1 / t116 ^ 2;
t41 = (-t117 * t95 * t186 + (t125 * t80 - t126 * t83) * t94) * t87 * t132;
t203 = t117 * t41;
t218 = rSges(3,1) * t116;
t38 = t41 ^ 2;
t29 = (-pkin(2) * t116 * t38 - t96 * t191) * t87;
t32 = t95 * t190 + t38;
t142 = -rSges(3,2) * t110 + t218;
t44 = (t142 * m(3) + t123) * t117 - t111 * t75;
t56 = g(1) * t80 + g(2) * t83;
t160 = t94 * t183;
t187 = t110 * t116;
t198 = t132 * t95;
t7 = ((pkin(2) * t93 * t203 - t111 * t186) * t41 * t198 + (-t111 * t41 * t187 + t117 * t160) * t96 * t183) * t87;
t221 = (-t44 * t7 - (t144 * t193 + t41 * t75 + t74) * t203 + (-t29 - t56) * t85 + (-t38 * t123 + (-t32 * t218 + (rSges(3,2) * t32 - t62 * t168) * t110) * m(3)) * t111) * t87;
t113 = sin(qJ(2,1));
t97 = t118 ^ 2;
t100 = t98 / t97;
t158 = t100 * t190;
t192 = t98 * t112;
t119 = cos(qJ(2,1));
t184 = t112 * t124;
t104 = legFrame(1,2);
t81 = sin(t104);
t84 = cos(t104);
t88 = 0.1e1 / t113;
t99 = 0.1e1 / t118 ^ 2;
t42 = (-t119 * t99 * t184 + (t125 * t81 - t126 * t84) * t98) * t88 * t132;
t202 = t119 * t42;
t217 = rSges(3,1) * t118;
t39 = t42 ^ 2;
t30 = (-pkin(2) * t118 * t39 - t100 * t191) * t88;
t33 = t99 * t190 + t39;
t141 = -rSges(3,2) * t112 + t217;
t45 = (t141 * m(3) + t123) * t119 - t113 * t75;
t57 = g(1) * t81 + g(2) * t84;
t159 = t98 * t183;
t185 = t112 * t118;
t196 = t132 * t99;
t8 = ((pkin(2) * t97 * t202 - t113 * t184) * t42 * t196 + (-t113 * t42 * t185 + t119 * t159) * t100 * t183) * t88;
t220 = (-t45 * t8 - (t144 * t192 + t42 * t75 + t74) * t202 + (-t30 - t57) * t85 + (-t39 * t123 + (-t33 * t217 + (rSges(3,2) * t33 - t63 * t158) * t112) * m(3)) * t113) * t88;
t106 = xDDP(2);
t213 = t106 * t86;
t212 = t106 * t87;
t211 = t106 * t88;
t107 = xDDP(1);
t210 = t107 * t86;
t209 = t107 * t87;
t208 = t107 * t88;
t207 = t108 * t86;
t206 = t110 * t87;
t205 = t112 * t88;
t201 = t132 * t90;
t199 = t132 * t94;
t197 = t132 * t98;
t182 = (t130 + t131);
t181 = t90 * t226;
t180 = t94 * t225;
t179 = t98 * t224;
t178 = t79 * t201;
t177 = t80 * t199;
t176 = t81 * t197;
t175 = t82 * t201;
t174 = t83 * t199;
t173 = t84 * t197;
t172 = t86 * t201;
t171 = t87 * t199;
t170 = t88 * t197;
t167 = t68 * t189;
t166 = t68 * t187;
t165 = t68 * t185;
t164 = t115 * t200;
t163 = t117 * t198;
t162 = t119 * t196;
t156 = t79 * t172;
t155 = t80 * t171;
t154 = t81 * t170;
t153 = t82 * t172;
t152 = t83 * t171;
t151 = t84 * t170;
t150 = t86 * t164;
t149 = t87 * t163;
t148 = t88 * t162;
t147 = t108 * t169;
t146 = t110 * t168;
t145 = t112 * t158;
t58 = g(1) * t82 - g(2) * t79;
t140 = t109 * t55 + t115 * t58;
t59 = g(1) * t83 - g(2) * t80;
t139 = t111 * t56 + t117 * t59;
t60 = g(1) * t84 - g(2) * t81;
t138 = t113 * t57 + t119 * t60;
t137 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (2 * rSges(3,3) ^ 2 + t182) * m(3) / 0.2e1;
t129 = 0.2e1 * qJ(3,1);
t128 = 0.2e1 * qJ(3,2);
t127 = 0.2e1 * qJ(3,3);
t122 = rSges(3,2) * g(3);
t105 = xDDP(3);
t73 = t182 * m(3) + Icges(3,3);
t54 = t113 * t81 + t119 * t84;
t53 = t113 * t84 - t119 * t81;
t52 = t111 * t80 + t117 * t83;
t51 = t111 * t83 - t117 * t80;
t50 = t109 * t79 + t115 * t82;
t49 = t109 * t82 - t115 * t79;
t48 = -t77 * t112 - t118 * t76;
t47 = -t77 * t110 - t116 * t76;
t46 = -t77 * t108 - t114 * t76;
t36 = cos(t129) * t229 - t78 * sin(t129) + t137;
t35 = cos(t128) * t229 - t78 * sin(t128) + t137;
t34 = cos(t127) * t229 - t78 * sin(t127) + t137;
t27 = (-t45 * t173 + t54 * t85) * t88;
t26 = (t45 * t176 + t53 * t85) * t88;
t25 = (-t44 * t174 + t52 * t85) * t87;
t24 = (t44 * t177 + t51 * t85) * t87;
t23 = (-t43 * t175 + t50 * t85) * t86;
t22 = (t43 * t178 + t49 * t85) * t86;
t21 = -t113 * t132 * t179 + (-t45 * t162 + t85 * t98) * t205;
t20 = -t111 * t132 * t180 + (-t44 * t163 + t85 * t94) * t206;
t19 = -t109 * t132 * t181 + (-t43 * t164 + t85 * t90) * t207;
t18 = (-t36 * t173 + t45 * t54) * t88;
t17 = (t36 * t176 + t45 * t53) * t88;
t16 = (-t35 * t174 + t44 * t52) * t87;
t15 = (t35 * t177 + t44 * t51) * t87;
t14 = (-t34 * t175 + t43 * t50) * t86;
t13 = (t34 * t178 + t43 * t49) * t86;
t12 = t48 * t197 + (-t36 * t162 + t45 * t98) * t205;
t11 = t47 * t199 + (-t35 * t163 + t44 * t94) * t206;
t10 = t46 * t201 + (-t34 * t164 + t43 * t90) * t207;
t3 = -t45 * t30 - t36 * t8 + t48 * t145 - 0.4e1 * ((t112 * t228 + t118 * t227) * t159 + (t165 / 0.2e1 + (t97 - 0.1e1 / 0.2e1) * t78) * t42) * t159 + (m(2) * (-rSges(2,1) * t57 + rSges(2,2) * t60) + (-t60 * rSges(3,3) - t141 * t57) * m(3)) * t119 + t113 * (m(2) * (rSges(2,1) * t60 + rSges(2,2) * t57) + (-t57 * rSges(3,3) + t141 * t60) * m(3));
t2 = -t44 * t29 - t35 * t7 + t47 * t146 - 0.4e1 * ((t110 * t228 + t116 * t227) * t160 + (t166 / 0.2e1 + (t93 - 0.1e1 / 0.2e1) * t78) * t41) * t160 + (m(2) * (-rSges(2,1) * t56 + rSges(2,2) * t59) + (-t59 * rSges(3,3) - t142 * t56) * m(3)) * t117 + t111 * (m(2) * (rSges(2,1) * t59 + rSges(2,2) * t56) + (-t56 * rSges(3,3) + t142 * t59) * m(3));
t1 = -t43 * t28 - t34 * t9 + t46 * t147 - 0.4e1 * ((t108 * t228 + t114 * t227) * t161 + (t167 / 0.2e1 + (t89 - 0.1e1 / 0.2e1) * t78) * t40) * t161 + (m(2) * (-rSges(2,1) * t55 + rSges(2,2) * t58) + (-t58 * rSges(3,3) - t143 * t55) * m(3)) * t115 + t109 * (m(2) * (rSges(2,1) * t58 + rSges(2,2) * t55) + (-t55 * rSges(3,3) + t143 * t58) * m(3));
t4 = [(-t18 * t173 + t27 * t54) * t208 + (t18 * t176 + t27 * t53) * t211 + t54 * t220 - t3 * t151 + (-t16 * t174 + t25 * t52) * t209 + (t16 * t177 + t25 * t51) * t212 + t52 * t221 - t2 * t152 + (-t14 * t175 + t23 * t50) * t210 + (t14 * t178 + t23 * t49) * t213 + t50 * t222 - t1 * t153 + (-g(1) + t107) * m(4) + ((-t48 * t151 - t54 * t224) * t197 + (-t18 * t162 + t27 * t98) * t205 + (-t47 * t152 - t52 * t225) * t199 + (-t16 * t163 + t25 * t94) * t206 + (-t46 * t153 - t50 * t226) * t201 + (-t14 * t164 + t23 * t90) * t207) * t105; (-t17 * t173 + t26 * t54) * t208 + (t17 * t176 + t26 * t53) * t211 + t53 * t220 + t3 * t154 + (-t15 * t174 + t24 * t52) * t209 + (t15 * t177 + t24 * t51) * t212 + t51 * t221 + t2 * t155 + (-t13 * t175 + t22 * t50) * t210 + (t13 * t178 + t22 * t49) * t213 + t49 * t222 + t1 * t156 + (-g(2) + t106) * m(4) + ((t48 * t154 - t53 * t224) * t197 + (-t17 * t162 + t26 * t98) * t205 + (t47 * t155 - t51 * t225) * t199 + (-t15 * t163 + t24 * t94) * t206 + (t46 * t156 - t49 * t226) * t201 + (-t13 * t164 + t22 * t90) * t207) * t105; (-t12 * t173 + t21 * t54) * t208 + (t12 * t176 + t21 * t53) * t211 + t192 * t220 - t112 * t3 * t148 + (-t48 * t8 + t73 * t145 + t39 * (t97 * t230 + Icges(3,4) + t165) + (t113 * t63 * t30 + t112 * t122 + t138 * t214 + (-t39 * rSges(3,2) - g(3) * t118 + t112 * t138) * rSges(3,1)) * m(3)) * t197 + (-t11 * t174 + t20 * t52) * t209 + (t11 * t177 + t20 * t51) * t212 + t193 * t221 - t110 * t2 * t149 + (-t47 * t7 + t73 * t146 + t38 * (t93 * t230 + Icges(3,4) + t166) + (t111 * t62 * t29 + t110 * t122 + t139 * t215 + (-t38 * rSges(3,2) - g(3) * t116 + t110 * t139) * rSges(3,1)) * m(3)) * t199 + (-t10 * t175 + t19 * t50) * t210 + (t10 * t178 + t19 * t49) * t213 + t194 * t222 - t108 * t1 * t150 + (-t46 * t9 + t73 * t147 + t37 * (t89 * t230 + Icges(3,4) + t167) + (t109 * t61 * t28 + t108 * t122 + t140 * t216 + (-t37 * rSges(3,2) - g(3) * t114 + t108 * t140) * rSges(3,1)) * m(3)) * t201 - m(4) * g(3) + ((t21 * t98 * t88 - t12 * t148 + (-t48 * t148 - t179) * t197) * t112 + (t20 * t94 * t87 - t11 * t149 + (-t47 * t149 - t180) * t199) * t110 + (t19 * t90 * t86 - t10 * t150 + (-t46 * t150 - t181) * t201) * t108 + m(4) + (t90 ^ 2 + t94 ^ 2 + t98 ^ 2) * t73 * t132 ^ 2) * t105;];
tauX  = t4;
