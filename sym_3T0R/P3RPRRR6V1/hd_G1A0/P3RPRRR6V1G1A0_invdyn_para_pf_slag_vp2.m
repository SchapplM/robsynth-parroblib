% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR6V1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:31:59
% EndTime: 2020-08-06 18:32:03
% DurationCPUTime: 4.03s
% Computational Cost: add. (15520->427), mult. (13597->601), div. (1527->13), fcn. (7575->98), ass. (0->298)
t196 = cos(pkin(7));
t167 = t196 * pkin(1);
t146 = t167 + pkin(2);
t349 = mrSges(3,1) * t146;
t195 = sin(pkin(7));
t308 = t195 * pkin(1);
t239 = pkin(5) + t308;
t348 = t239 * mrSges(3,1);
t347 = t239 * mrSges(3,2);
t217 = 2 * qJ(3,3);
t187 = sin(t217);
t158 = pkin(3) * t187;
t204 = sin(qJ(3,3));
t248 = pkin(7) + qJ(3,3);
t251 = -pkin(7) + qJ(3,3);
t343 = 0.2e1 * pkin(2);
t70 = 0.1e1 / (t204 * t343 + t158 + (sin(t248) + sin(t251)) * pkin(1));
t346 = t70 / 0.2e1;
t218 = 2 * qJ(3,2);
t188 = sin(t218);
t159 = pkin(3) * t188;
t205 = sin(qJ(3,2));
t249 = pkin(7) + qJ(3,2);
t252 = -pkin(7) + qJ(3,2);
t71 = 0.1e1 / (t205 * t343 + t159 + (sin(t249) + sin(t252)) * pkin(1));
t345 = t71 / 0.2e1;
t219 = 2 * qJ(3,1);
t189 = sin(t219);
t160 = pkin(3) * t189;
t206 = sin(qJ(3,1));
t250 = pkin(7) + qJ(3,1);
t253 = -pkin(7) + qJ(3,1);
t72 = 0.1e1 / (t206 * t343 + t160 + (sin(t250) + sin(t253)) * pkin(1));
t344 = t72 / 0.2e1;
t333 = -0.2e1 * pkin(2);
t186 = qJ(1,1) + pkin(7);
t199 = legFrame(1,3);
t139 = t199 + t186;
t120 = cos(t139);
t130 = t219 + t139;
t133 = -(2 * qJ(3,1)) + t139;
t157 = qJ(1,1) + t199;
t144 = qJ(3,1) + t157;
t145 = -qJ(3,1) + t157;
t256 = 0.2e1 * pkin(1);
t131 = qJ(3,1) + t139;
t132 = -qJ(3,1) + t139;
t260 = cos(t131) + cos(t132);
t263 = sin(t131) + sin(t132);
t214 = -pkin(6) - pkin(5);
t321 = -0.2e1 * t214;
t54 = t260 * t343 + (cos(t145) + cos(t144)) * t256 + t263 * t321 + (cos(t133) + cos(t130) + 0.2e1 * t120) * pkin(3);
t340 = t54 * t344;
t185 = qJ(1,2) + pkin(7);
t198 = legFrame(2,3);
t138 = t198 + t185;
t119 = cos(t138);
t126 = t218 + t138;
t129 = -(2 * qJ(3,2)) + t138;
t156 = qJ(1,2) + t198;
t142 = qJ(3,2) + t156;
t143 = -qJ(3,2) + t156;
t163 = qJ(1,2) + t249;
t127 = t198 + t163;
t164 = qJ(1,2) - t252;
t128 = t198 + t164;
t261 = cos(t127) + cos(t128);
t264 = sin(t127) + sin(t128);
t53 = t261 * t343 + (cos(t143) + cos(t142)) * t256 + t264 * t321 + (cos(t129) + cos(t126) + 0.2e1 * t119) * pkin(3);
t339 = t53 * t345;
t184 = qJ(1,3) + pkin(7);
t197 = legFrame(3,3);
t137 = t197 + t184;
t118 = cos(t137);
t122 = t217 + t137;
t125 = -(2 * qJ(3,3)) + t137;
t155 = qJ(1,3) + t197;
t140 = qJ(3,3) + t155;
t141 = -qJ(3,3) + t155;
t123 = qJ(3,3) + t137;
t124 = -qJ(3,3) + t137;
t262 = cos(t123) + cos(t124);
t265 = sin(t123) + sin(t124);
t52 = t262 * t343 + (cos(t141) + cos(t140)) * t256 + t265 * t321 + (cos(t125) + cos(t122) + 0.2e1 * t118) * pkin(3);
t338 = t52 * t346;
t117 = sin(t139);
t320 = 0.2e1 * t214;
t51 = t263 * t343 + (sin(t145) + sin(t144)) * t256 + t260 * t320 + (sin(t133) + sin(t130) + 0.2e1 * t117) * pkin(3);
t337 = t51 * t344;
t116 = sin(t138);
t50 = t264 * t343 + (sin(t143) + sin(t142)) * t256 + t261 * t320 + (sin(t129) + sin(t126) + 0.2e1 * t116) * pkin(3);
t336 = t50 * t345;
t115 = sin(t137);
t49 = t265 * t343 + (sin(t141) + sin(t140)) * t256 + t262 * t320 + (sin(t125) + sin(t122) + 0.2e1 * t115) * pkin(3);
t335 = t49 * t346;
t334 = -0.2e1 * pkin(1);
t212 = xDP(2);
t213 = xDP(1);
t221 = 0.1e1 / pkin(3);
t280 = t221 * t70;
t58 = t118 * t321 + sin(t155) * t334 + t115 * t333 - t265 * pkin(3);
t61 = t115 * t320 + cos(t155) * t334 + t118 * t333 - t262 * pkin(3);
t46 = (t212 * t58 + t213 * t61) * t280;
t43 = t46 ^ 2;
t279 = t221 * t71;
t59 = t119 * t321 + sin(t156) * t334 + t116 * t333 - t264 * pkin(3);
t62 = t116 * t320 + cos(t156) * t334 + t119 * t333 - t261 * pkin(3);
t47 = (t212 * t59 + t213 * t62) * t279;
t44 = t47 ^ 2;
t278 = t221 * t72;
t60 = t120 * t321 + sin(t157) * t334 + t117 * t333 - t263 * pkin(3);
t63 = t117 * t320 + cos(t157) * t334 + t120 * t333 - t260 * pkin(3);
t48 = (t212 * t60 + t213 * t63) * t278;
t45 = t48 ^ 2;
t207 = cos(qJ(3,3));
t294 = t207 * pkin(3) + pkin(2);
t93 = t167 + t294;
t89 = 0.1e1 / t93;
t64 = (-t115 * t213 + t118 * t212) * t89;
t331 = t64 ^ 2;
t208 = cos(qJ(3,2));
t293 = t208 * pkin(3) + pkin(2);
t94 = t167 + t293;
t90 = 0.1e1 / t94;
t65 = (-t116 * t213 + t119 * t212) * t90;
t330 = t65 ^ 2;
t209 = cos(qJ(3,1));
t292 = t209 * pkin(3) + pkin(2);
t95 = t167 + t292;
t91 = 0.1e1 / t95;
t66 = (-t117 * t213 + t120 * t212) * t91;
t329 = t66 ^ 2;
t220 = pkin(3) ^ 2;
t222 = pkin(2) ^ 2;
t223 = pkin(1) ^ 2;
t233 = t214 ^ 2 + t308 * t321 + t222 + t223;
t254 = 0.2e1 * t167;
t328 = t254 * t333 - t220 - 0.2e1 * t233;
t191 = t207 ^ 2;
t327 = -0.2e1 * t191;
t192 = t208 ^ 2;
t326 = -0.2e1 * t192;
t193 = t209 ^ 2;
t325 = -0.2e1 * t193;
t324 = 0.2e1 * t207;
t323 = 0.2e1 * t208;
t322 = 0.2e1 * t209;
t216 = m(2) + m(3);
t319 = mrSges(3,2) * pkin(2);
t168 = sin(t197);
t171 = cos(t197);
t83 = -t168 * g(1) + t171 * g(2);
t318 = mrSges(3,1) * t83;
t169 = sin(t198);
t172 = cos(t198);
t84 = -t169 * g(1) + t172 * g(2);
t317 = mrSges(3,1) * t84;
t170 = sin(t199);
t173 = cos(t199);
t85 = -t170 * g(1) + t173 * g(2);
t316 = mrSges(3,1) * t85;
t315 = mrSges(3,2) * t83;
t314 = mrSges(3,2) * t84;
t313 = mrSges(3,2) * t85;
t86 = t171 * g(1) + t168 * g(2);
t312 = mrSges(3,2) * t86;
t87 = t172 * g(1) + t169 * g(2);
t311 = mrSges(3,2) * t87;
t88 = t173 * g(1) + t170 * g(2);
t310 = mrSges(3,2) * t88;
t309 = g(3) * t216;
t301 = t58 * t70;
t300 = t59 * t71;
t299 = t60 * t72;
t298 = t61 * t70;
t297 = t62 * t71;
t296 = t63 * t72;
t291 = Ifges(3,3) * t221;
t290 = t115 * t89;
t289 = t116 * t90;
t288 = t117 * t91;
t287 = t118 * t89;
t286 = t119 * t90;
t285 = t120 * t91;
t284 = t206 * mrSges(3,2);
t97 = -t209 * mrSges(3,1) + t284;
t283 = t97 / 0.2e1;
t176 = t204 * mrSges(3,2);
t98 = t207 * mrSges(3,1) - t176;
t282 = t98 / 0.2e1;
t177 = t205 * mrSges(3,2);
t99 = t208 * mrSges(3,1) - t177;
t281 = t99 / 0.2e1;
t277 = t221 * t97;
t276 = t221 * t98;
t275 = t221 * t99;
t148 = m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3);
t274 = t148 * t195;
t200 = Ifges(3,1) - Ifges(3,2);
t273 = t200 * t204;
t272 = t200 * t205;
t271 = t200 * t206;
t270 = t204 * t146;
t269 = t205 * t146;
t268 = t206 * t146;
t267 = t216 / 0.2e1;
t266 = pkin(3) * t343;
t257 = -0.2e1 * t319;
t247 = t46 * t204 * pkin(3);
t246 = t47 * t205 * pkin(3);
t245 = t48 * t206 * pkin(3);
t31 = (t49 * t267 + t58 * t276) * t70;
t32 = (t50 * t267 + t59 * t275) * t71;
t33 = (t51 * t267 - t60 * t277) * t72;
t244 = t33 + t32 + t31;
t34 = (t52 * t267 + t61 * t276) * t70;
t35 = (t53 * t267 + t62 * t275) * t71;
t36 = (t54 * t267 - t63 * t277) * t72;
t243 = t36 + t35 + t34;
t225 = Ifges(3,5) - t348;
t96 = Ifges(3,6) - t347;
t67 = t204 * t225 + t96 * t207;
t242 = t67 * t280;
t68 = t205 * t225 + t96 * t208;
t241 = t68 * t279;
t69 = t206 * t225 + t96 * t209;
t240 = t69 * t278;
t179 = m(3) * pkin(2) + mrSges(2,1);
t238 = (t327 + 0.1e1) * Ifges(3,4);
t237 = (t326 + 0.1e1) * Ifges(3,4);
t236 = (t325 + 0.1e1) * Ifges(3,4);
t235 = -t214 + t308;
t100 = mrSges(3,1) * t204 + mrSges(3,2) * t207;
t135 = pkin(6) + t239;
t74 = 0.1e1 / (t158 + 0.2e1 * t270);
t13 = (t331 * t207 * t328 + (-0.2e1 * t43 * t93 + ((-cos(t217) * t93 + t146 * t327 - t146) * t64 + (t204 * t324 + t187) * t46 * t135) * t64) * pkin(3)) * t74;
t19 = (-t135 * t247 + (t191 * t220 + t207 * t266 + t294 * t254 + t233) * t64) / (t158 / 0.2e1 + t270) * t64 * t221 + 0.2e1 * (-t64 * t135 * t204 + t207 * t46 * t93) * t74 * t46;
t232 = -t43 * t100 - t216 * t13 - t98 * t19;
t101 = mrSges(3,1) * t205 + mrSges(3,2) * t208;
t75 = 0.1e1 / (t159 + 0.2e1 * t269);
t14 = (t330 * t208 * t328 + (-0.2e1 * t44 * t94 + ((-cos(t218) * t94 + t146 * t326 - t146) * t65 + (t205 * t323 + t188) * t47 * t135) * t65) * pkin(3)) * t75;
t20 = (-t135 * t246 + (t192 * t220 + t208 * t266 + t293 * t254 + t233) * t65) / (t159 / 0.2e1 + t269) * t65 * t221 + 0.2e1 * (-t65 * t135 * t205 + t208 * t47 * t94) * t75 * t47;
t231 = -t44 * t101 - t216 * t14 - t99 * t20;
t102 = mrSges(3,1) * t206 + mrSges(3,2) * t209;
t76 = 0.1e1 / (t160 + 0.2e1 * t268);
t15 = (t329 * t209 * t328 + (-0.2e1 * t45 * t95 + ((-cos(t219) * t95 + t146 * t325 - t146) * t66 + (t206 * t322 + t189) * t48 * t135) * t66) * pkin(3)) * t76;
t21 = (-t135 * t245 + (t193 * t220 + t209 * t266 + t292 * t254 + t233) * t66) / (t160 / 0.2e1 + t268) * t66 * t221 + 0.2e1 * (-t66 * t135 * t206 + t209 * t48 * t95) * t76 * t48;
t230 = -t45 * t102 - t216 * t15 + t97 * t21;
t228 = Ifges(3,6) / 0.2e1 - t347 / 0.2e1;
t227 = -Ifges(3,5) / 0.2e1 + t348 / 0.2e1;
t226 = mrSges(3,2) * t167 + t319;
t224 = (pkin(5) ^ 2 + t222) * m(3) + 0.2e1 * mrSges(3,3) * pkin(5) + t216 * t223 + Ifges(3,1) + Ifges(1,3) + Ifges(2,3);
t203 = xDDP(1);
t202 = xDDP(2);
t201 = xDDP(3);
t166 = qJ(1,1) - t253;
t165 = qJ(1,1) + t250;
t162 = qJ(1,3) - t251;
t161 = qJ(1,3) + t248;
t154 = cos(t186);
t153 = cos(t185);
t152 = cos(t184);
t151 = sin(t186);
t150 = sin(t185);
t149 = sin(t184);
t136 = t216 * pkin(1) + mrSges(1,1);
t121 = mrSges(3,2) * t146;
t79 = mrSges(3,1) * t88;
t78 = mrSges(3,1) * t87;
t77 = mrSges(3,1) * t86;
t57 = -t200 * t193 + (Ifges(3,4) * t206 + t349) * t322 + t206 * t257 + (-(-t179 + t284) * t196 + t274) * t256 + t224;
t56 = -t200 * t192 + (Ifges(3,4) * t205 + t349) * t323 + t205 * t257 + (-(t177 - t179) * t196 + t274) * t256 + t224;
t55 = -t200 * t191 + (Ifges(3,4) * t204 + t349) * t324 + t204 * t257 + (-(t176 - t179) * t196 + t274) * t256 + t224;
t42 = t60 * t240 + t57 * t285;
t41 = t59 * t241 + t56 * t286;
t40 = t58 * t242 + t55 * t287;
t39 = t63 * t240 - t57 * t288;
t38 = t62 * t241 - t56 * t289;
t37 = t61 * t242 - t55 * t290;
t30 = t69 * t285 + (-t51 * t283 + t60 * t291) * t72;
t29 = t68 * t286 + (t50 * t281 + t59 * t291) * t71;
t28 = t67 * t287 + (t49 * t282 + t58 * t291) * t70;
t27 = -t69 * t288 + (-t54 * t283 + t63 * t291) * t72;
t26 = -t68 * t289 + (t53 * t281 + t62 * t291) * t71;
t25 = -t67 * t290 + (t52 * t282 + t61 * t291) * t70;
t24 = (t235 * t65 - 0.2e1 * t246) * t90 * t65;
t23 = (t235 * t64 - 0.2e1 * t247) * t89 * t64;
t22 = (t235 * t66 - 0.2e1 * t245) * t91 * t66;
t9 = t230 - t309;
t8 = t231 - t309;
t7 = t232 - t309;
t6 = -t69 * t22 - Ifges(3,3) * t21 + t329 * ((t121 - t271) * t209 + mrSges(3,1) * t268 + t236) + (t15 + g(3)) * t97 + (t85 * t151 + t88 * t154) * t102;
t5 = -t68 * t24 - Ifges(3,3) * t20 + t330 * ((t121 - t272) * t208 + mrSges(3,1) * t269 + t237) + (-t14 - g(3)) * t99 + (t84 * t150 + t87 * t153) * t101;
t4 = -t67 * t23 - Ifges(3,3) * t19 + t331 * ((t121 - t273) * t207 + mrSges(3,1) * t270 + t238) + (-t13 - g(3)) * t98 + (t83 * t149 + t86 * t152) * t100;
t3 = -t57 * t22 - t69 * t21 - 0.2e1 * ((t228 * t206 + t227 * t209) * t48 + ((t226 - t271) * t209 + t206 * t349 + t236) * t66) * t48 + (-t310 - t316) * cos(t166) / 0.2e1 + (t79 - t313) * sin(t166) / 0.2e1 + (t310 - t316) * cos(t165) / 0.2e1 + (t79 + t313) * sin(t165) / 0.2e1 + (mrSges(1,2) * t88 - t85 * t136) * cos(qJ(1,1)) + (t85 * mrSges(1,2) + t136 * t88) * sin(qJ(1,1)) + (-t148 * t88 - t85 * t179) * t154 + (-t85 * t148 + t179 * t88) * t151;
t2 = -t56 * t24 - t68 * t20 - 0.2e1 * ((t228 * t205 + t227 * t208) * t47 + ((t226 - t272) * t208 + t205 * t349 + t237) * t65) * t47 + (-t311 - t317) * cos(t164) / 0.2e1 + (t78 - t314) * sin(t164) / 0.2e1 + (t311 - t317) * cos(t163) / 0.2e1 + (t78 + t314) * sin(t163) / 0.2e1 + (mrSges(1,2) * t87 - t84 * t136) * cos(qJ(1,2)) + (t84 * mrSges(1,2) + t136 * t87) * sin(qJ(1,2)) + (-t148 * t87 - t84 * t179) * t153 + (-t84 * t148 + t179 * t87) * t150;
t1 = -t55 * t23 - t67 * t19 - 0.2e1 * ((t228 * t204 + t227 * t207) * t46 + ((t226 - t273) * t207 + t204 * t349 + t238) * t64) * t46 + (-t312 - t318) * cos(t162) / 0.2e1 + (t77 - t315) * sin(t162) / 0.2e1 + (t312 - t318) * cos(t161) / 0.2e1 + (t77 + t315) * sin(t161) / 0.2e1 + (mrSges(1,2) * t86 - t83 * t136) * cos(qJ(1,3)) + (t83 * mrSges(1,2) + t136 * t86) * sin(qJ(1,3)) + (-t148 * t86 - t83 * t179) * t152 + (-t83 * t148 + t179 * t86) * t149;
t10 = [-t1 * t290 - t2 * t289 - t3 * t288 - g(1) * m(4) + (t6 * t296 + t5 * t297 + t4 * t298) * t221 + t9 * t340 + t8 * t339 + t7 * t338 + (t37 * t287 + t38 * t286 + t39 * t285 + (t25 * t301 + t26 * t300 + t27 * t299) * t221 + t36 * t337 + t35 * t336 + t34 * t335) * t202 + t243 * t201 + (-t37 * t290 - t38 * t289 - t39 * t288 + m(4) + (t25 * t298 + t26 * t297 + t27 * t296) * t221 + t36 * t340 + t35 * t339 + t34 * t338) * t203; t1 * t287 + t2 * t286 + t3 * t285 - g(2) * m(4) + (t6 * t299 + t5 * t300 + t4 * t301) * t221 + t9 * t337 + t8 * t336 + t7 * t335 + (-t40 * t290 - t41 * t289 - t42 * t288 + (t28 * t298 + t29 * t297 + t30 * t296) * t221 + t33 * t340 + t32 * t339 + t31 * t338) * t203 + t244 * t201 + (t40 * t287 + t41 * t286 + t42 * t285 + m(4) + (t28 * t301 + t29 * t300 + t30 * t299) * t221 + t33 * t337 + t32 * t336 + t31 * t335) * t202; t244 * t202 + t243 * t203 + t230 + t231 + t232 + (-g(3) + t201) * (m(4) + 0.3e1 * t216);];
tauX  = t10;
