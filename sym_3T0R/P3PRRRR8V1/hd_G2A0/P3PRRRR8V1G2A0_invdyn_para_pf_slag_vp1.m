% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:02:29
% EndTime: 2020-08-06 17:02:35
% DurationCPUTime: 6.82s
% Computational Cost: add. (18822->401), mult. (58296->796), div. (7479->9), fcn. (61890->28), ass. (0->327)
t181 = cos(qJ(3,3));
t182 = cos(qJ(2,3));
t176 = sin(qJ(2,3));
t279 = t176 * t181;
t128 = pkin(2) * t279 - pkin(5) * t182;
t166 = sin(pkin(3));
t168 = cos(pkin(3));
t175 = sin(qJ(3,3));
t110 = pkin(2) * t175 * t168 + t128 * t166;
t373 = 0.1e1 / t110;
t326 = t373 / t181;
t171 = legFrame(1,2);
t154 = sin(t171);
t157 = cos(t171);
t124 = g(1) * t154 + g(2) * t157;
t127 = g(1) * t157 - g(2) * t154;
t165 = sin(pkin(6));
t151 = g(3) * t165;
t167 = cos(pkin(6));
t241 = t127 * t167 - t151;
t384 = t124 * t166 + t241 * t168;
t170 = legFrame(2,2);
t153 = sin(t170);
t156 = cos(t170);
t123 = g(1) * t153 + g(2) * t156;
t126 = g(1) * t156 - g(2) * t153;
t243 = t126 * t167 - t151;
t383 = t123 * t166 + t243 * t168;
t169 = legFrame(3,2);
t152 = sin(t169);
t155 = cos(t169);
t122 = g(1) * t152 + g(2) * t155;
t125 = g(1) * t155 - g(2) * t152;
t245 = t125 * t167 - t151;
t382 = t122 * t166 + t245 * t168;
t179 = sin(qJ(3,1));
t185 = cos(qJ(3,1));
t378 = -rSges(3,1) * t185 + rSges(3,2) * t179;
t177 = sin(qJ(3,2));
t183 = cos(qJ(3,2));
t377 = -rSges(3,1) * t183 + rSges(3,2) * t177;
t376 = -rSges(3,1) * t181 + rSges(3,2) * t175;
t190 = xDP(3);
t191 = xDP(2);
t192 = xDP(1);
t237 = t152 * t191 - t155 * t192;
t289 = t168 * t176;
t213 = t166 * t181 + t175 * t289;
t280 = t175 * t182;
t88 = t213 * t165 - t167 * t280;
t91 = -t165 * t280 - t213 * t167;
t64 = (t190 * t91 + t237 * t88) * t326;
t61 = t64 ^ 2;
t236 = t153 * t191 - t156 * t192;
t162 = 0.1e1 / t183;
t178 = sin(qJ(2,2));
t277 = t178 * t183;
t251 = t166 * t277;
t288 = t168 * t177;
t184 = cos(qJ(2,2));
t291 = t166 * t184;
t375 = 0.1e1 / (-pkin(5) * t291 + (t251 + t288) * pkin(2));
t325 = t375 * t162;
t287 = t168 * t178;
t212 = t166 * t183 + t177 * t287;
t278 = t177 * t184;
t89 = t212 * t165 - t167 * t278;
t92 = -t165 * t278 - t212 * t167;
t65 = (t190 * t92 + t236 * t89) * t325;
t62 = t65 ^ 2;
t235 = t154 * t191 - t157 * t192;
t164 = 0.1e1 / t185;
t180 = sin(qJ(2,1));
t275 = t180 * t185;
t250 = t166 * t275;
t286 = t168 * t179;
t186 = cos(qJ(2,1));
t290 = t166 * t186;
t374 = 0.1e1 / (-pkin(5) * t290 + (t250 + t286) * pkin(2));
t324 = t374 * t164;
t285 = t168 * t180;
t211 = t166 * t185 + t179 * t285;
t276 = t179 * t186;
t90 = t211 * t165 - t167 * t276;
t93 = -t165 * t276 - t211 * t167;
t66 = (t190 * t93 + t235 * t90) * t324;
t63 = t66 ^ 2;
t149 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t372 = 0.2e1 * t149;
t189 = m(2) * rSges(2,1);
t371 = m(3) * rSges(3,3);
t343 = rSges(3,2) * t181;
t135 = rSges(3,1) * t175 + t343;
t94 = -t166 * t176 * t135 - t168 * t376;
t370 = m(3) * t94;
t342 = rSges(3,2) * t183;
t136 = rSges(3,1) * t177 + t342;
t95 = -t166 * t178 * t136 - t168 * t377;
t369 = m(3) * t95;
t341 = rSges(3,2) * t185;
t137 = rSges(3,1) * t179 + t341;
t96 = -t166 * t180 * t137 - t168 * t378;
t368 = m(3) * t96;
t200 = 0.1e1 / pkin(2);
t284 = t168 * t182;
t354 = pkin(2) * t181;
t82 = -(-t165 * t176 + t167 * t284) * t354 - pkin(5) * (t165 * t182 + t167 * t289);
t85 = (t165 * t284 + t167 * t176) * t354 + (t165 * t289 - t167 * t182) * pkin(5);
t58 = (t190 * t82 + t237 * t85) * t200 * t326;
t367 = pkin(2) * t58;
t283 = t168 * t184;
t353 = pkin(2) * t183;
t83 = -(-t165 * t178 + t167 * t283) * t353 - pkin(5) * (t165 * t184 + t167 * t287);
t86 = (t165 * t283 + t167 * t178) * t353 + (t165 * t287 - t167 * t184) * pkin(5);
t59 = (t190 * t83 + t236 * t86) * t200 * t325;
t366 = pkin(2) * t59;
t282 = t168 * t186;
t352 = pkin(2) * t185;
t84 = -(-t165 * t180 + t167 * t282) * t352 - pkin(5) * (t165 * t186 + t167 * t285);
t87 = (t165 * t282 + t167 * t180) * t352 + (t165 * t285 - t167 * t186) * pkin(5);
t60 = (t190 * t84 + t235 * t87) * t200 * t324;
t365 = pkin(2) * t60;
t364 = pkin(5) * t64;
t363 = pkin(5) * t65;
t362 = pkin(5) * t66;
t196 = rSges(3,2) ^ 2;
t197 = rSges(3,1) ^ 2;
t139 = (-t196 + t197) * m(3) - Icges(3,1) + Icges(3,2);
t361 = t139 / 0.2e1;
t147 = rSges(3,2) * t371 - Icges(3,6);
t360 = -t147 / 0.4e1;
t148 = rSges(3,1) * t371 - Icges(3,5);
t359 = t148 / 0.4e1;
t358 = m(3) * t200;
t159 = t181 ^ 2;
t357 = pkin(2) * t159;
t161 = t183 ^ 2;
t356 = pkin(2) * t161;
t163 = t185 ^ 2;
t355 = pkin(2) * t163;
t351 = g(3) * t167;
t198 = pkin(5) ^ 2;
t199 = pkin(2) ^ 2;
t335 = t175 * t58;
t267 = pkin(2) * t335;
t350 = (-pkin(5) * t267 + (t159 * t199 + t198) * t64) * t64;
t340 = t375 * t65;
t339 = t374 * t66;
t119 = -m(3) * t376 + t189;
t146 = m(2) * rSges(2,2) - t371;
t97 = t119 * t182 - t146 * t176;
t338 = t166 * t97;
t120 = -m(3) * t377 + t189;
t98 = t120 * t184 - t146 * t178;
t337 = t166 * t98;
t121 = -m(3) * t378 + t189;
t99 = t121 * t186 - t146 * t180;
t336 = t166 * t99;
t334 = t177 * t59;
t333 = t179 * t60;
t332 = t200 * t82;
t331 = t200 * t83;
t330 = t200 * t84;
t329 = t200 * t85;
t328 = t200 * t86;
t327 = t200 * t87;
t158 = m(1) + m(2) + m(3);
t323 = t373 * t158;
t129 = pkin(2) * t277 - pkin(5) * t184;
t111 = pkin(2) * t288 + t129 * t166;
t108 = 0.1e1 / t111;
t322 = t108 * t158;
t130 = pkin(2) * t275 - pkin(5) * t186;
t112 = pkin(2) * t286 + t130 * t166;
t109 = 0.1e1 / t112;
t321 = t109 * t158;
t116 = -t147 * t181 - t148 * t175;
t320 = t116 * t200;
t117 = -t147 * t183 - t148 * t177;
t319 = t117 * t200;
t118 = -t147 * t185 - t148 * t179;
t318 = t118 * t200;
t316 = t122 * t168;
t314 = t123 * t168;
t312 = t124 * t168;
t308 = t139 * t175;
t307 = t139 * t177;
t306 = t139 * t179;
t274 = t196 + t197;
t144 = t274 * m(3) + Icges(3,3);
t305 = t144 * t200;
t304 = t146 * t165;
t303 = rSges(3,2) * t151 * t166;
t173 = xDDP(2);
t302 = t152 * t173;
t301 = t153 * t173;
t300 = t154 * t173;
t174 = xDDP(1);
t299 = t155 * t174;
t298 = t156 * t174;
t297 = t157 * t174;
t296 = t166 * t167;
t295 = t166 * t175;
t294 = t166 * t177;
t293 = t166 * t179;
t292 = t166 * t182;
t281 = t168 * t200;
t273 = t373 * t370;
t272 = t108 * t369;
t271 = t109 * t368;
t270 = t94 * t358;
t269 = t95 * t358;
t268 = t96 * t358;
t266 = pkin(2) * t334;
t265 = pkin(2) * t333;
t264 = t175 * t364;
t263 = t177 * t363;
t262 = t179 * t362;
t261 = t373 * t338;
t260 = t108 * t337;
t259 = t109 * t336;
t258 = t152 * t326;
t257 = t155 * t326;
t256 = t153 * t325;
t255 = t156 * t325;
t254 = t154 * t324;
t253 = t157 * t324;
t252 = t166 * t279;
t246 = t125 * t165 + t351;
t244 = t126 * t165 + t351;
t242 = t127 * t165 + t351;
t22 = t264 - t367;
t10 = (((t168 * t58 + t64 * t292) * t357 - (t267 - t364) * t252 + t168 * t22) * t64 - (-t58 * t292 + (-t159 * t168 + t175 * t252 + t168) * t64) * t367) * t326;
t13 = (t281 * t350 + (-t58 * t128 * t295 + t168 * (t58 * t357 - t264)) * t58) * t326;
t134 = t146 * t351;
t16 = (t22 * t367 - t350) * t373;
t193 = 0.2e1 * qJ(3,3);
t204 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(3,3) ^ 2 + t274) * m(3) / 0.2e1;
t73 = cos(t193) * t361 - t149 * sin(t193) + t204;
t4 = -t16 * t338 - t73 * t10 - t116 * t13 - 0.4e1 * t58 * ((t64 * t308 / 0.2e1 + t58 * t359) * t181 + t335 * t360 + (t159 - 0.1e1 / 0.2e1) * t64 * t149) + (-t119 * t382 + t125 * t304 + t134) * t182 - t176 * (-t246 * t119 - t146 * t382);
t203 = t176 * t382 + t246 * t182;
t7 = -t16 * t370 - t116 * t10 - t144 * t13 + (t159 * t372 + t181 * t308 - t149) * t61 + m(3) * (((t245 * t166 - t316) * rSges(3,1) + t203 * rSges(3,2)) * t181 + t175 * (t303 + (-t125 * t296 + t316) * rSges(3,2) + t203 * rSges(3,1)));
t234 = t7 * t329 + t88 * t4;
t23 = t263 - t366;
t11 = (((t168 * t59 + t65 * t291) * t356 - (t266 - t363) * t251 + t168 * t23) * t340 + (t59 * t291 + (t161 * t168 - t177 * t251 - t168) * t65) * t375 * t366) * t162;
t20 = -pkin(5) * t266 + (t161 * t199 + t198) * t65;
t14 = (t20 * t281 * t340 + (-t59 * t129 * t294 + t168 * (t59 * t356 - t263)) * t108 * t59) * t162;
t17 = (t20 * t65 - t23 * t366) * t375;
t194 = 0.2e1 * qJ(3,2);
t74 = cos(t194) * t361 - t149 * sin(t194) + t204;
t5 = t17 * t337 - t74 * t11 - t117 * t14 - 0.4e1 * t59 * ((t65 * t307 / 0.2e1 + t59 * t359) * t183 + t334 * t360 + (t161 - 0.1e1 / 0.2e1) * t65 * t149) + (-t120 * t383 + t126 * t304 + t134) * t184 - t178 * (-t244 * t120 - t146 * t383);
t202 = t178 * t383 + t244 * t184;
t8 = t17 * t369 - t117 * t11 - t144 * t14 + (t161 * t372 + t183 * t307 - t149) * t62 + m(3) * (((t243 * t166 - t314) * rSges(3,1) + t202 * rSges(3,2)) * t183 + t177 * (t303 + (-t126 * t296 + t314) * rSges(3,2) + t202 * rSges(3,1)));
t233 = t8 * t328 + t89 * t5;
t24 = t262 - t365;
t12 = (((t168 * t60 + t66 * t290) * t355 - (t265 - t362) * t250 + t168 * t24) * t339 + (t60 * t290 + (t163 * t168 - t179 * t250 - t168) * t66) * t374 * t365) * t164;
t21 = -pkin(5) * t265 + (t163 * t199 + t198) * t66;
t15 = (t21 * t281 * t339 + (-t60 * t130 * t293 + t168 * (t60 * t355 - t262)) * t109 * t60) * t164;
t18 = (t21 * t66 - t24 * t365) * t374;
t195 = 0.2e1 * qJ(3,1);
t75 = cos(t195) * t361 - t149 * sin(t195) + t204;
t6 = t18 * t336 - t75 * t12 - t118 * t15 - 0.4e1 * t60 * ((t66 * t306 / 0.2e1 + t60 * t359) * t185 + t333 * t360 + (t163 - 0.1e1 / 0.2e1) * t66 * t149) + (-t121 * t384 + t127 * t304 + t134) * t186 - t180 * (-t242 * t121 - t146 * t384);
t201 = t180 * t384 + t242 * t186;
t9 = t18 * t368 - t118 * t12 - t144 * t15 + (t163 * t372 + t185 * t306 - t149) * t63 + m(3) * (((t241 * t166 - t312) * rSges(3,1) + t201 * rSges(3,2)) * t185 + t179 * (t303 + (-t127 * t296 + t312) * rSges(3,2) + t201 * rSges(3,1)));
t232 = t9 * t327 + t90 * t6;
t222 = t85 * t320 + t73 * t88;
t131 = pkin(5) * t176 + t182 * t354;
t219 = pkin(2) * t295 - t128 * t168;
t76 = -t165 * t131 + t219 * t167;
t67 = t110 * t155 + t152 * t76;
t28 = t222 * t258 + t67 * t261;
t216 = t116 * t88 + t85 * t305;
t40 = t216 * t258 + t67 * t273;
t231 = t28 * t88 + t40 * t329;
t68 = t110 * t152 - t155 * t76;
t29 = -t222 * t257 + t68 * t261;
t41 = -t216 * t257 + t68 * t273;
t230 = t29 * t88 + t41 * t329;
t221 = t86 * t319 + t74 * t89;
t132 = pkin(5) * t178 + t184 * t353;
t218 = pkin(2) * t294 - t129 * t168;
t77 = -t165 * t132 + t218 * t167;
t69 = t111 * t156 + t153 * t77;
t30 = t221 * t256 + t69 * t260;
t215 = t117 * t89 + t86 * t305;
t42 = t215 * t256 + t69 * t272;
t229 = t30 * t89 + t42 * t328;
t70 = t111 * t153 - t156 * t77;
t31 = -t221 * t255 + t70 * t260;
t43 = -t215 * t255 + t70 * t272;
t228 = t31 * t89 + t43 * t328;
t220 = t87 * t318 + t75 * t90;
t133 = pkin(5) * t180 + t186 * t352;
t217 = pkin(2) * t293 - t130 * t168;
t78 = -t165 * t133 + t217 * t167;
t71 = t112 * t157 + t154 * t78;
t32 = t220 * t254 + t71 * t259;
t214 = t118 * t90 + t87 * t305;
t44 = t214 * t254 + t71 * t271;
t227 = t32 * t90 + t44 * t327;
t72 = t112 * t154 - t157 * t78;
t33 = -t220 * t253 + t72 * t259;
t45 = -t214 * t253 + t72 * t271;
t226 = t45 * t327 + t33 * t90;
t210 = t85 * t270 + t88 * t338;
t209 = t86 * t269 + t89 * t337;
t208 = t87 * t268 + t90 * t336;
t172 = xDDP(3);
t79 = t167 * t131 + t219 * t165;
t207 = t172 * t79 + t173 * t67 + t174 * t68;
t80 = t167 * t132 + t218 * t165;
t206 = t172 * t80 + t173 * t69 + t174 * t70;
t81 = t167 * t133 + t217 * t165;
t205 = t172 * t81 + t173 * t71 + t174 * t72;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = t81 * t271 + (t118 * t93 + t84 * t305) * t324;
t53 = t80 * t272 + (t117 * t92 + t83 * t305) * t325;
t52 = t79 * t273 + (t116 * t91 + t82 * t305) * t326;
t48 = t81 * t259 + (t84 * t318 + t75 * t93) * t324;
t47 = t80 * t260 + (t83 * t319 + t74 * t92) * t325;
t46 = t79 * t261 + (t82 * t320 + t73 * t91) * t326;
t3 = (-t12 * t99 + (-t146 * t186 - t180 * t189) * t63) * t166 + (t18 - t124) * t158 + (-t96 * t15 + (-0.2e1 * (rSges(3,1) * t333 + t60 * t341) * t66 * t186 + t378 * t180 * (t63 + t57)) * t166 - t57 * t168 * t137) * m(3);
t2 = (-t11 * t98 + (-t146 * t184 - t178 * t189) * t62) * t166 + (t17 - t123) * t158 + (-t95 * t14 + (-0.2e1 * (rSges(3,1) * t334 + t59 * t342) * t65 * t184 + t377 * t178 * (t62 + t56)) * t166 - t56 * t168 * t136) * m(3);
t1 = (-t10 * t97 + (-t146 * t182 - t176 * t189) * t61) * t166 + (-t16 - t122) * t158 + (-t94 * t13 + (-0.2e1 * (rSges(3,1) * t335 + t58 * t343) * t64 * t182 + t376 * t176 * (t61 + t55)) * t166 - t55 * t168 * t135) * m(3);
t19 = [(-g(1) + t174) * m(4) + (t72 * t3 + t205 * (-t208 * t253 + t72 * t321)) * t109 + (t70 * t2 + t206 * (-t209 * t255 + t70 * t322)) * t108 + (t68 * t1 + t207 * (-t210 * t257 + t68 * t323)) * t373 + ((t33 * t93 + t45 * t330) * t172 + t226 * t300 + (-t226 * t174 - t232) * t157) * t324 + ((t31 * t92 + t43 * t331) * t172 + t228 * t301 + (-t228 * t174 - t233) * t156) * t325 + ((t29 * t91 + t41 * t332) * t172 + t230 * t302 + (-t230 * t174 - t234) * t155) * t326; (-g(2) + t173) * m(4) + (t71 * t3 + t205 * (t208 * t254 + t71 * t321)) * t109 + (t69 * t2 + t206 * (t209 * t256 + t69 * t322)) * t108 + (t67 * t1 + t207 * (t210 * t258 + t67 * t323)) * t373 + ((t32 * t93 + t44 * t330) * t172 - t227 * t297 + (t227 * t173 + t232) * t154) * t324 + ((t30 * t92 + t42 * t331) * t172 - t229 * t298 + (t229 * t173 + t233) * t153) * t325 + ((t28 * t91 + t40 * t332) * t172 - t231 * t299 + (t231 * t173 + t234) * t152) * t326; (-g(3) + t172) * m(4) + (t81 * t3 + t205 * (t81 * t321 + (t84 * t268 + t93 * t336) * t324)) * t109 + (t80 * t2 + t206 * (t80 * t322 + (t83 * t269 + t92 * t337) * t325)) * t108 + (t79 * t1 + t207 * (t79 * t323 + (t82 * t270 + t91 * t338) * t326)) * t373 + ((t54 * t330 + t48 * t93) * t172 + t93 * t6 + t9 * t330 + (-t297 + t300) * (t54 * t327 + t48 * t90)) * t324 + ((t53 * t331 + t47 * t92) * t172 + t92 * t5 + t8 * t331 + (-t298 + t301) * (t53 * t328 + t47 * t89)) * t325 + ((t52 * t332 + t46 * t91) * t172 + t91 * t4 + t7 * t332 + (-t299 + t302) * (t52 * t329 + t46 * t88)) * t326;];
tauX  = t19;
