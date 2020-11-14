% Calculate inertia matrix for parallel robot
% P3PRRR1G1A0
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
% koppelP [3x3]
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:50
% EndTime: 2020-03-09 21:14:50
% DurationCPUTime: 0.22s
% Computational Cost: add. (1153->75), mult. (735->154), div. (180->5), fcn. (792->24), ass. (0->91)
t261 = pkin(7) + qJ(2,3);
t252 = qJ(3,3) + t261;
t240 = sin(t252);
t243 = cos(t252);
t246 = sin(t261);
t249 = cos(t261);
t309 = 0.1e1 / (t240 * t249 - t246 * t243);
t262 = pkin(7) + qJ(2,2);
t253 = qJ(3,2) + t262;
t241 = sin(t253);
t244 = cos(t253);
t247 = sin(t262);
t250 = cos(t262);
t308 = 0.1e1 / (t241 * t250 - t247 * t244);
t263 = pkin(7) + qJ(2,1);
t254 = qJ(3,1) + t263;
t242 = sin(t254);
t245 = cos(t254);
t248 = sin(t263);
t251 = cos(t263);
t307 = 0.1e1 / (t242 * t251 - t248 * t245);
t264 = legFrame(3,3);
t255 = sin(t264);
t258 = cos(t264);
t228 = t258 * t240 + t255 * t243;
t213 = pkin(2) * (t246 * t258 + t255 * t249) + t228 * pkin(3);
t306 = t213 * t309;
t265 = legFrame(2,3);
t256 = sin(t265);
t259 = cos(t265);
t230 = t259 * t241 + t256 * t244;
t214 = pkin(2) * (t247 * t259 + t256 * t250) + t230 * pkin(3);
t305 = t214 * t308;
t266 = legFrame(1,3);
t257 = sin(t266);
t260 = cos(t266);
t232 = t260 * t242 + t257 * t245;
t215 = pkin(2) * (t248 * t260 + t257 * t251) + t232 * pkin(3);
t304 = t215 * t307;
t229 = -t240 * t255 + t258 * t243;
t216 = -pkin(2) * (t246 * t255 - t258 * t249) + t229 * pkin(3);
t303 = t216 * t309;
t231 = -t241 * t256 + t259 * t244;
t217 = -pkin(2) * (t247 * t256 - t259 * t250) + t231 * pkin(3);
t302 = t217 * t308;
t233 = -t242 * t257 + t260 * t245;
t218 = -pkin(2) * (t248 * t257 - t260 * t251) + t233 * pkin(3);
t301 = t218 * t307;
t300 = t309 * t228;
t299 = t309 * t229;
t278 = (mrSges(3,1) * cos(qJ(3,3)) - mrSges(3,2) * sin(qJ(3,3))) * pkin(2);
t285 = m(3) * pkin(2) ^ 2 + Ifges(2,3) + Ifges(3,3);
t298 = t309 * (0.2e1 * t278 + t285);
t237 = Ifges(3,3) + t278;
t297 = t309 * t237;
t296 = t308 * t230;
t295 = t308 * t231;
t277 = (mrSges(3,1) * cos(qJ(3,2)) - mrSges(3,2) * sin(qJ(3,2))) * pkin(2);
t294 = t308 * (0.2e1 * t277 + t285);
t238 = Ifges(3,3) + t277;
t293 = t308 * t238;
t292 = t307 * t232;
t291 = t307 * t233;
t276 = (mrSges(3,1) * cos(qJ(3,1)) - mrSges(3,2) * sin(qJ(3,1))) * pkin(2);
t290 = t307 * (0.2e1 * t276 + t285);
t239 = Ifges(3,3) + t276;
t289 = t307 * t239;
t274 = 0.1e1 / pkin(3);
t288 = t309 * t274;
t287 = t308 * t274;
t286 = t307 * t274;
t284 = Ifges(3,3) * t288;
t283 = Ifges(3,3) * t287;
t282 = Ifges(3,3) * t286;
t281 = t237 * t288;
t280 = t238 * t287;
t279 = t239 * t286;
t275 = 0.1e1 / pkin(2);
t212 = (-t218 * t282 + t233 * t289) * t275;
t211 = (-t215 * t282 + t232 * t289) * t275;
t210 = (-t217 * t283 + t231 * t293) * t275;
t209 = (-t214 * t283 + t230 * t293) * t275;
t208 = (-t216 * t284 + t229 * t297) * t275;
t207 = (-t213 * t284 + t228 * t297) * t275;
t206 = (-t218 * t279 + t233 * t290) * t275;
t205 = (-t215 * t279 + t232 * t290) * t275;
t204 = (-t217 * t280 + t231 * t294) * t275;
t203 = (-t214 * t280 + t230 * t294) * t275;
t202 = (-t216 * t281 + t229 * t298) * t275;
t201 = (-t213 * t281 + t228 * t298) * t275;
t1 = [m(4) + (t202 * t299 + t204 * t295 + t206 * t291 + (-t208 * t303 - t210 * t302 - t212 * t301) * t274) * t275, (t202 * t300 + t204 * t296 + t206 * t292 + (-t208 * t306 - t210 * t305 - t212 * t304) * t274) * t275, 0; (t201 * t299 + t203 * t295 + t205 * t291 + (-t207 * t303 - t209 * t302 - t211 * t301) * t274) * t275, m(4) + (t201 * t300 + t203 * t296 + t205 * t292 + (-t207 * t306 - t209 * t305 - t211 * t304) * t274) * t275, 0; 0, 0, (3 * m(1)) + (3 * m(2)) + 0.3e1 * m(3) + m(4);];
MX  = t1;
