% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR2G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:04
% EndTime: 2020-03-09 21:20:04
% DurationCPUTime: 0.39s
% Computational Cost: add. (1768->87), mult. (1290->174), div. (1212->5), fcn. (1260->18), ass. (0->107)
t299 = cos(qJ(3,3));
t347 = t299 * pkin(1);
t300 = cos(qJ(3,2));
t346 = t300 * pkin(1);
t301 = cos(qJ(3,1));
t345 = t301 * pkin(1);
t303 = xDP(1);
t305 = 0.1e1 / pkin(2);
t307 = 0.1e1 / pkin(1);
t321 = t305 * t307;
t315 = t303 * t321;
t287 = -legFrame(3,2) + qJ(2,3);
t284 = qJ(3,3) + t287;
t278 = sin(t284);
t263 = pkin(2) * t278 + pkin(1) * sin(t287);
t296 = sin(qJ(3,3));
t290 = 0.1e1 / t296;
t341 = t263 * t290;
t257 = t315 * t341;
t302 = xDP(2);
t316 = t302 * t321;
t281 = cos(t284);
t266 = -pkin(2) * t281 - pkin(1) * cos(t287);
t338 = t266 * t290;
t260 = t316 * t338;
t248 = t257 + t260;
t322 = t303 * t307;
t323 = t302 * t307;
t329 = t281 * t290;
t332 = t278 * t290;
t254 = -t322 * t332 + t323 * t329;
t245 = t254 + t248;
t344 = t245 * t248;
t288 = -legFrame(2,2) + qJ(2,2);
t285 = qJ(3,2) + t288;
t279 = sin(t285);
t264 = pkin(2) * t279 + pkin(1) * sin(t288);
t297 = sin(qJ(3,2));
t291 = 0.1e1 / t297;
t340 = t264 * t291;
t258 = t315 * t340;
t282 = cos(t285);
t267 = -pkin(2) * t282 - pkin(1) * cos(t288);
t337 = t267 * t291;
t261 = t316 * t337;
t249 = t258 + t261;
t328 = t282 * t291;
t331 = t279 * t291;
t255 = -t322 * t331 + t323 * t328;
t246 = t255 + t249;
t343 = t246 * t249;
t289 = -legFrame(1,2) + qJ(2,1);
t286 = qJ(3,1) + t289;
t280 = sin(t286);
t265 = pkin(2) * t280 + pkin(1) * sin(t289);
t298 = sin(qJ(3,1));
t292 = 0.1e1 / t298;
t339 = t265 * t292;
t259 = t315 * t339;
t283 = cos(t286);
t268 = -pkin(2) * t283 - pkin(1) * cos(t289);
t336 = t268 * t292;
t262 = t316 * t336;
t250 = t259 + t262;
t327 = t283 * t292;
t330 = t280 * t292;
t256 = -t322 * t330 + t323 * t327;
t247 = t256 + t250;
t342 = t247 * t250;
t335 = (t296 * mrSges(3,1) + t299 * mrSges(3,2)) * t290;
t334 = (t297 * mrSges(3,1) + t300 * mrSges(3,2)) * t291;
t333 = (t298 * mrSges(3,1) + t301 * mrSges(3,2)) * t292;
t326 = t290 * t307;
t325 = t291 * t307;
t324 = t292 * t307;
t320 = 0.2e1 * pkin(1) * pkin(2);
t319 = t254 ^ 2 * t335;
t318 = t255 ^ 2 * t334;
t317 = t256 ^ 2 * t333;
t242 = t257 / 0.2e1 + t260 / 0.2e1 + t254;
t314 = t242 * t248 * t335;
t243 = t258 / 0.2e1 + t261 / 0.2e1 + t255;
t313 = t243 * t249 * t334;
t244 = t259 / 0.2e1 + t262 / 0.2e1 + t256;
t312 = t244 * t250 * t333;
t306 = pkin(1) ^ 2;
t311 = -m(3) * t306 - Ifges(2,3) - Ifges(3,3);
t310 = (-mrSges(3,1) * t299 + mrSges(3,2) * t296) * pkin(1);
t309 = (-mrSges(3,1) * t300 + mrSges(3,2) * t297) * pkin(1);
t308 = (-mrSges(3,1) * t301 + mrSges(3,2) * t298) * pkin(1);
t304 = pkin(2) ^ 2;
t274 = -Ifges(3,3) + t308;
t273 = -Ifges(3,3) + t309;
t272 = -Ifges(3,3) + t310;
t241 = ((-t247 * pkin(2) - t256 * t345) * t256 - pkin(2) * t342) * t324;
t240 = ((-pkin(2) * t246 - t255 * t346) * t255 - pkin(2) * t343) * t325;
t239 = ((-pkin(2) * t245 - t254 * t347) * t254 - pkin(2) * t344) * t326;
t238 = ((t243 * t300 * t320 + t246 * t304 + t306 * t255) * t305 * t255 + (pkin(2) + t346) * t343) * t325;
t237 = ((t242 * t299 * t320 + t245 * t304 + t306 * t254) * t305 * t254 + (pkin(2) + t347) * t344) * t326;
t236 = ((t244 * t301 * t320 + t247 * t304 + t306 * t256) * t305 * t256 + (pkin(2) + t345) * t342) * t324;
t235 = -Ifges(3,3) * t236 + t274 * t241;
t234 = -Ifges(3,3) * t238 + t273 * t240;
t233 = -Ifges(3,3) * t237 + t272 * t239;
t232 = (0.2e1 * t308 + t311) * t241 + t274 * t236;
t231 = (0.2e1 * t309 + t311) * t240 + t273 * t238;
t230 = (0.2e1 * t310 + t311) * t239 + t272 * t237;
t1 = [0.2e1 * t278 * t314 + 0.2e1 * t279 * t313 + 0.2e1 * t280 * t312 + (-t230 * t332 - t231 * t331 - t232 * t330) * t307 + (t263 * t319 + t264 * t318 + t265 * t317 + (t233 * t341 + t234 * t340 + t235 * t339) * t307) * t305; -0.2e1 * t281 * t314 - 0.2e1 * t282 * t313 - 0.2e1 * t283 * t312 + (t230 * t329 + t231 * t328 + t232 * t327) * t307 + (t266 * t319 + t267 * t318 + t268 * t317 + (t233 * t338 + t234 * t337 + t235 * t336) * t307) * t305; 0;];
taucX  = t1;
