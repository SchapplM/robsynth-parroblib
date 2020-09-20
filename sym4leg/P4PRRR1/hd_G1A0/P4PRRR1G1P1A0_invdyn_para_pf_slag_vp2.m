% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% tauX [4x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:12:24
% EndTime: 2020-03-02 20:12:30
% DurationCPUTime: 6.56s
% Computational Cost: add. (21885->313), mult. (14553->589), div. (3184->6), fcn. (15454->34), ass. (0->265)
t229 = 0.1e1 / pkin(2);
t190 = pkin(7) + qJ(2,4);
t169 = qJ(3,4) + t190;
t155 = sin(t169);
t156 = cos(t169);
t167 = sin(t190);
t168 = cos(t190);
t319 = 0.1e1 / (t155 * t168 - t156 * t167);
t323 = t229 * t319;
t193 = pkin(7) + qJ(2,1);
t178 = qJ(3,1) + t193;
t159 = sin(t178);
t162 = cos(t178);
t172 = sin(t193);
t175 = cos(t193);
t316 = 0.1e1 / (t159 * t175 - t162 * t172);
t322 = t316 * t229;
t192 = pkin(7) + qJ(2,2);
t177 = qJ(3,2) + t192;
t158 = sin(t177);
t161 = cos(t177);
t171 = sin(t192);
t174 = cos(t192);
t317 = 0.1e1 / (t158 * t174 - t161 * t171);
t321 = t317 * t229;
t191 = pkin(7) + qJ(2,3);
t176 = qJ(3,3) + t191;
t157 = sin(t176);
t160 = cos(t176);
t170 = sin(t191);
t173 = cos(t191);
t318 = 0.1e1 / (t157 * t173 - t160 * t170);
t320 = t318 * t229;
t315 = mrSges(3,1) * pkin(2);
t314 = mrSges(3,2) * pkin(2);
t200 = sin(qJ(3,4));
t201 = cos(qJ(3,4));
t313 = pkin(2) * (mrSges(3,1) * t200 + mrSges(3,2) * t201);
t206 = sin(qJ(3,3));
t209 = cos(qJ(3,3));
t312 = pkin(2) * (mrSges(3,1) * t206 + mrSges(3,2) * t209);
t207 = sin(qJ(3,2));
t210 = cos(qJ(3,2));
t311 = pkin(2) * (mrSges(3,1) * t207 + mrSges(3,2) * t210);
t208 = sin(qJ(3,1));
t211 = cos(qJ(3,1));
t310 = pkin(2) * (mrSges(3,1) * t208 + mrSges(3,2) * t211);
t163 = t200 * t314;
t228 = pkin(2) ^ 2;
t246 = m(3) * t228 + Ifges(2,3) + Ifges(3,3);
t271 = t201 * t315;
t125 = -0.2e1 * t163 + t246 + 0.2e1 * t271;
t139 = Ifges(3,3) - t163 + t271;
t215 = xP(4);
t188 = sin(t215);
t189 = cos(t215);
t218 = koppelP(4,2);
t222 = koppelP(4,1);
t131 = t188 * t222 + t189 * t218;
t135 = -t188 * t218 + t189 * t222;
t227 = 0.1e1 / pkin(3);
t272 = t227 * t229;
t258 = t319 * t272;
t196 = legFrame(4,3);
t179 = sin(t196);
t183 = cos(t196);
t109 = t155 * t183 + t156 * t179;
t73 = pkin(2) * (t167 * t183 + t168 * t179) + t109 * pkin(3);
t110 = -t155 * t179 + t156 * t183;
t74 = -pkin(2) * (t167 * t179 - t168 * t183) + t110 * pkin(3);
t37 = (-t131 * t74 + t135 * t73) * t258;
t61 = (t109 * t135 - t110 * t131) * t323;
t309 = (t125 * t61 - t139 * t37) * t319;
t212 = xDP(4);
t214 = xDP(1);
t121 = -t131 * t212 + t214;
t294 = t74 * t121;
t213 = xDP(2);
t117 = t135 * t212 + t213;
t295 = t73 * t117;
t233 = (t294 + t295) * t258;
t49 = (t109 * t117 + t110 * t121) * t323;
t19 = -t233 + t49;
t308 = t19 * t233;
t255 = t318 * t272;
t219 = koppelP(3,2);
t223 = koppelP(3,1);
t132 = t188 * t223 + t189 * t219;
t122 = -t132 * t212 + t214;
t197 = legFrame(3,3);
t180 = sin(t197);
t184 = cos(t197);
t112 = -t157 * t180 + t160 * t184;
t78 = -pkin(2) * (t170 * t180 - t173 * t184) + t112 * pkin(3);
t290 = t78 * t122;
t136 = -t188 * t219 + t189 * t223;
t118 = t136 * t212 + t213;
t111 = t157 * t184 + t160 * t180;
t75 = pkin(2) * (t170 * t184 + t173 * t180) + t111 * pkin(3);
t293 = t75 * t118;
t232 = (t290 + t293) * t255;
t52 = (t111 * t118 + t112 * t122) * t320;
t26 = -t232 + t52;
t307 = t26 * t232;
t253 = t317 * t272;
t220 = koppelP(2,2);
t224 = koppelP(2,1);
t133 = t188 * t224 + t189 * t220;
t123 = -t133 * t212 + t214;
t198 = legFrame(2,3);
t181 = sin(t198);
t185 = cos(t198);
t114 = -t158 * t181 + t161 * t185;
t79 = -pkin(2) * (t171 * t181 - t174 * t185) + t114 * pkin(3);
t289 = t79 * t123;
t137 = -t188 * t220 + t189 * t224;
t119 = t137 * t212 + t213;
t113 = t158 * t185 + t161 * t181;
t76 = pkin(2) * (t171 * t185 + t174 * t181) + t113 * pkin(3);
t292 = t76 * t119;
t231 = (t289 + t292) * t253;
t53 = (t113 * t119 + t114 * t123) * t321;
t27 = -t231 + t53;
t306 = t27 * t231;
t251 = t316 * t272;
t221 = koppelP(1,2);
t225 = koppelP(1,1);
t134 = t188 * t225 + t189 * t221;
t124 = -t134 * t212 + t214;
t199 = legFrame(1,3);
t182 = sin(t199);
t186 = cos(t199);
t116 = -t159 * t182 + t162 * t186;
t80 = -pkin(2) * (t172 * t182 - t175 * t186) + t116 * pkin(3);
t288 = t80 * t124;
t138 = -t188 * t221 + t189 * t225;
t120 = t138 * t212 + t213;
t115 = t159 * t186 + t162 * t182;
t77 = pkin(2) * (t172 * t186 + t175 * t182) + t115 * pkin(3);
t291 = t77 * t120;
t230 = (t288 + t291) * t251;
t54 = (t115 * t120 + t116 * t124) * t322;
t28 = -t230 + t54;
t305 = t28 * t230;
t164 = t206 * t314;
t270 = t209 * t315;
t126 = -0.2e1 * t164 + t246 + 0.2e1 * t270;
t148 = Ifges(3,3) - t164 + t270;
t38 = (-t132 * t78 + t136 * t75) * t255;
t62 = (t111 * t136 - t112 * t132) * t320;
t303 = t318 * (t126 * t62 - t148 * t38);
t165 = t207 * t314;
t269 = t210 * t315;
t127 = -0.2e1 * t165 + t246 + 0.2e1 * t269;
t149 = Ifges(3,3) - t165 + t269;
t39 = (-t133 * t79 + t137 * t76) * t253;
t63 = (t113 * t137 - t114 * t133) * t321;
t302 = t317 * (t127 * t63 - t149 * t39);
t166 = t208 * t314;
t268 = t211 * t315;
t128 = -0.2e1 * t166 + t246 + 0.2e1 * t268;
t150 = Ifges(3,3) - t166 + t268;
t40 = (-t134 * t80 + t138 * t77) * t251;
t64 = (t115 * t138 - t116 * t134) * t322;
t301 = t316 * (t128 * t64 - t150 * t40);
t300 = t125 * t319;
t299 = t139 * t319;
t298 = t227 * t319;
t287 = t318 * t126;
t286 = t318 * t148;
t284 = t317 * t127;
t283 = t317 * t149;
t281 = t316 * t128;
t280 = t316 * t150;
t278 = t318 * t227;
t276 = t317 * t227;
t274 = t316 * t227;
t267 = t73 * t298;
t266 = t74 * t298;
t265 = t75 * t278;
t264 = t78 * t278;
t263 = t76 * t276;
t262 = t79 * t276;
t261 = t77 * t274;
t260 = t80 * t274;
t259 = t139 * t298;
t257 = 0.2e1 * pkin(2) * pkin(3);
t256 = t148 * t278;
t254 = t149 * t276;
t252 = t150 * t274;
t140 = -g(1) * t179 + g(2) * t183;
t144 = g(1) * t183 + g(2) * t179;
t250 = -t156 * (mrSges(3,1) * t140 - mrSges(3,2) * t144) + (mrSges(3,1) * t144 + mrSges(3,2) * t140) * t155;
t141 = -g(1) * t180 + g(2) * t184;
t145 = g(1) * t184 + g(2) * t180;
t249 = -t160 * (mrSges(3,1) * t141 - mrSges(3,2) * t145) + (mrSges(3,1) * t145 + mrSges(3,2) * t141) * t157;
t142 = -g(1) * t181 + g(2) * t185;
t146 = g(1) * t185 + g(2) * t181;
t248 = -t161 * (mrSges(3,1) * t142 - mrSges(3,2) * t146) + (mrSges(3,1) * t146 + mrSges(3,2) * t142) * t158;
t143 = -g(1) * t182 + g(2) * t186;
t147 = g(1) * t186 + g(2) * t182;
t247 = -t162 * (mrSges(3,1) * t143 - mrSges(3,2) * t147) + (mrSges(3,1) * t147 + mrSges(3,2) * t143) * t159;
t245 = t155 * t167 + t156 * t168;
t244 = t157 * t170 + t160 * t173;
t243 = t158 * t171 + t161 * t174;
t242 = t159 * t172 + t162 * t175;
t195 = t212 ^ 2;
t202 = xDDP(4);
t204 = xDDP(2);
t81 = -t131 * t195 + t135 * t202 + t204;
t205 = xDDP(1);
t85 = -t131 * t202 - t135 * t195 + t205;
t241 = t229 * (t73 * t81 + t74 * t85);
t82 = -t132 * t195 + t136 * t202 + t204;
t86 = -t132 * t202 - t136 * t195 + t205;
t240 = t229 * (t75 * t82 + t78 * t86);
t83 = -t133 * t195 + t137 * t202 + t204;
t87 = -t133 * t202 - t137 * t195 + t205;
t239 = t229 * (t76 * t83 + t79 * t87);
t84 = -t134 * t195 + t138 * t202 + t204;
t88 = -t134 * t202 - t138 * t195 + t205;
t238 = t229 * (t77 * t84 + t80 * t88);
t237 = t245 * pkin(2);
t236 = t244 * pkin(2);
t235 = t243 * pkin(2);
t234 = t242 * pkin(2);
t226 = pkin(3) ^ 2;
t217 = mrSges(4,1);
t216 = mrSges(4,2);
t187 = m(3) * pkin(2) + mrSges(2,1);
t130 = -t188 * t216 + t189 * t217;
t129 = t188 * t217 + t189 * t216;
t48 = (t116 * t281 - t252 * t80) * t229;
t47 = (t115 * t281 - t252 * t77) * t229;
t46 = (t114 * t284 - t254 * t79) * t229;
t45 = (t113 * t284 - t254 * t76) * t229;
t44 = (t112 * t287 - t256 * t78) * t229;
t43 = (t111 * t287 - t256 * t75) * t229;
t42 = (t110 * t300 - t259 * t74) * t229;
t41 = (t109 * t300 - t259 * t73) * t229;
t32 = -Ifges(3,3) * t40 + t150 * t64;
t31 = -Ifges(3,3) * t39 + t149 * t63;
t30 = -Ifges(3,3) * t38 + t148 * t62;
t29 = -Ifges(3,3) * t37 + t139 * t61;
t25 = -(t288 / 0.2e1 + t291 / 0.2e1) * t251 + t54;
t24 = -(t289 / 0.2e1 + t292 / 0.2e1) * t253 + t53;
t23 = -(t290 / 0.2e1 + t293 / 0.2e1) * t255 + t52;
t18 = -(t294 / 0.2e1 + t295 / 0.2e1) * t258 + t49;
t16 = (-pkin(3) * t305 + (pkin(3) * t28 + t234 * t54) * t54) * t322;
t15 = (-pkin(3) * t306 + (pkin(3) * t27 + t235 * t53) * t53) * t321;
t14 = (-pkin(3) * t307 + (pkin(3) * t26 + t236 * t52) * t52) * t320;
t13 = (-pkin(3) * t308 + (pkin(3) * t19 + t237 * t49) * t49) * t323;
t12 = ((-t242 * t25 * t257 - t226 * t28 - t228 * t54) * t227 * t54 + (pkin(3) + t234) * t305) * t322;
t11 = ((-t24 * t243 * t257 - t226 * t27 - t228 * t53) * t227 * t53 + (pkin(3) + t235) * t306) * t321;
t10 = ((-t23 * t244 * t257 - t226 * t26 - t228 * t52) * t227 * t52 + (pkin(3) + t236) * t307) * t320;
t9 = ((-t18 * t245 * t257 - t19 * t226 - t228 * t49) * t227 * t49 + (pkin(3) + t237) * t308) * t323;
t8 = t310 * t54 ^ 2 + Ifges(3,3) * t12 + t150 * t16 + t247;
t7 = t311 * t53 ^ 2 + Ifges(3,3) * t11 + t149 * t15 + t248;
t6 = t312 * t52 ^ 2 + Ifges(3,3) * t10 + t14 * t148 + t249;
t5 = t313 * t49 ^ 2 + Ifges(3,3) * t9 + t13 * t139 + t250;
t4 = t128 * t16 + t150 * t12 + 0.2e1 * t230 * t25 * t310 + (mrSges(2,2) * t147 - t143 * t187) * t175 + (mrSges(2,2) * t143 + t147 * t187) * t172 + t247;
t3 = t127 * t15 + t149 * t11 + 0.2e1 * t231 * t24 * t311 + (mrSges(2,2) * t146 - t142 * t187) * t174 + (mrSges(2,2) * t142 + t146 * t187) * t171 + t248;
t2 = t126 * t14 + t148 * t10 + 0.2e1 * t232 * t23 * t312 + (mrSges(2,2) * t145 - t141 * t187) * t173 + (mrSges(2,2) * t141 + t145 * t187) * t170 + t249;
t1 = t125 * t13 + t139 * t9 + 0.2e1 * t233 * t18 * t313 + (mrSges(2,2) * t144 - t140 * t187) * t168 + (mrSges(2,2) * t140 + t144 * t187) * t167 + t250;
t17 = [-t129 * t202 - t195 * t130 + (t205 - g(1)) * m(4) + ((t109 * t42 * t81 + (t42 * t85 + t1) * t110) * t319 + (t115 * t48 * t84 + (t48 * t88 + t4) * t116) * t316 + (t113 * t46 * t83 + (t46 * t87 + t3) * t114) * t317 + (t111 * t44 * t82 + (t44 * t86 + t2) * t112) * t318 + (-(t5 * t74 + (-Ifges(3,3) * t266 + t110 * t299) * t241) * t319 - (t8 * t80 + (-Ifges(3,3) * t260 + t116 * t280) * t238) * t316 - (t7 * t79 + (-Ifges(3,3) * t262 + t114 * t283) * t239) * t317 - (t6 * t78 + (-Ifges(3,3) * t264 + t112 * t286) * t240) * t318) * t227) * t229; -t195 * t129 + t130 * t202 + (t204 - g(2)) * m(4) + ((t110 * t41 * t85 + (t41 * t81 + t1) * t109) * t319 + (t116 * t47 * t88 + (t47 * t84 + t4) * t115) * t316 + (t114 * t45 * t87 + (t45 * t83 + t3) * t113) * t317 + (t112 * t43 * t86 + (t43 * t82 + t2) * t111) * t318 + (-(t5 * t73 + (-Ifges(3,3) * t267 + t109 * t299) * t241) * t319 - (t77 * t8 + (-Ifges(3,3) * t261 + t115 * t280) * t238) * t316 - (t7 * t76 + (-Ifges(3,3) * t263 + t113 * t283) * t239) * t317 - (t6 * t75 + (-Ifges(3,3) * t265 + t111 * t286) * t240) * t318) * t227) * t229; (-g(3) + xDDP(3)) * ((4 * m(1)) + (4 * m(2)) + 0.4e1 * m(3) + m(4)); t64 * t4 - t40 * t8 + t63 * t3 - t39 * t7 + t62 * t2 - t38 * t6 + t61 * t1 - t37 * t5 - t129 * t205 + t130 * t204 + Ifges(4,3) * t202 - (-g(1) * t217 - g(2) * t216) * t188 + t189 * (g(1) * t216 - g(2) * t217) + ((t116 * t301 - t260 * t32) * t88 + (t115 * t301 - t261 * t32) * t84 + (t114 * t302 - t262 * t31) * t87 + (t113 * t302 - t263 * t31) * t83 + (t112 * t303 - t264 * t30) * t86 + (t111 * t303 - t265 * t30) * t82 + (t110 * t309 - t266 * t29) * t85 + (t109 * t309 - t29 * t267) * t81) * t229;];
tauX  = t17;
