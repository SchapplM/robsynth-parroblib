% Calculate Gravitation load for parallel robot
% P3RPRR1G1P1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G1P1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:12
% EndTime: 2020-03-09 21:23:12
% DurationCPUTime: 0.24s
% Computational Cost: add. (380->89), mult. (497->124), div. (18->4), fcn. (288->48), ass. (0->57)
t523 = m(1) * rSges(1,2);
t522 = m(2) * rSges(2,2);
t521 = m(3) / pkin(3);
t509 = legFrame(2,3);
t500 = sin(t509);
t503 = cos(t509);
t465 = -t500 * g(1) + t503 * g(2);
t468 = t503 * g(1) + t500 * g(2);
t513 = pkin(7) + qJ(3,2);
t496 = qJ(1,2) + t513;
t455 = -cos(t496) * (rSges(3,1) * t465 - rSges(3,2) * t468) + sin(t496) * (rSges(3,1) * t468 + rSges(3,2) * t465);
t471 = 0.1e1 / (pkin(1) * sin(t513) + sin(qJ(3,2)) * pkin(2));
t520 = t455 * t471;
t508 = legFrame(3,3);
t499 = sin(t508);
t502 = cos(t508);
t464 = -t499 * g(1) + t502 * g(2);
t467 = t502 * g(1) + t499 * g(2);
t512 = pkin(7) + qJ(3,3);
t497 = qJ(1,3) + t512;
t456 = -cos(t497) * (rSges(3,1) * t464 - rSges(3,2) * t467) + sin(t497) * (rSges(3,1) * t467 + rSges(3,2) * t464);
t470 = 0.1e1 / (pkin(1) * sin(t512) + sin(qJ(3,3)) * pkin(2));
t519 = t456 * t470;
t510 = legFrame(1,3);
t501 = sin(t510);
t504 = cos(t510);
t466 = -t501 * g(1) + t504 * g(2);
t469 = t504 * g(1) + t501 * g(2);
t514 = pkin(7) + qJ(3,1);
t498 = qJ(1,1) + t514;
t457 = -cos(t498) * (rSges(3,1) * t466 - rSges(3,2) * t469) + sin(t498) * (rSges(3,1) * t469 + rSges(3,2) * t466);
t472 = 0.1e1 / (pkin(1) * sin(t514) + sin(qJ(3,1)) * pkin(2));
t518 = t457 * t472;
t479 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t492 = m(2) * rSges(2,1) + m(3) * pkin(2);
t505 = qJ(1,3) + pkin(7);
t517 = t470 * ((-t464 * t492 + t467 * t522) * cos(t505) + (t464 * t522 + t492 * t467) * sin(t505) + (-t464 * t479 + t467 * t523) * cos(qJ(1,3)) + sin(qJ(1,3)) * (t464 * t523 + t479 * t467) + t456 * m(3));
t506 = qJ(1,2) + pkin(7);
t516 = t471 * ((-t465 * t492 + t468 * t522) * cos(t506) + (t465 * t522 + t492 * t468) * sin(t506) + (-t465 * t479 + t468 * t523) * cos(qJ(1,2)) + sin(qJ(1,2)) * (t465 * t523 + t479 * t468) + t455 * m(3));
t507 = qJ(1,1) + pkin(7);
t515 = t472 * ((-t466 * t492 + t469 * t522) * cos(t507) + (t466 * t522 + t492 * t469) * sin(t507) + (-t466 * t479 + t469 * t523) * cos(qJ(1,1)) + sin(qJ(1,1)) * (t466 * t523 + t479 * t469) + t457 * m(3));
t485 = t510 + t507;
t483 = t508 + t505;
t495 = t510 + qJ(1,1);
t494 = t509 + qJ(1,2);
t493 = t508 + qJ(1,3);
t484 = t509 + t506;
t482 = qJ(3,1) + t485;
t481 = t509 + t496;
t480 = qJ(3,3) + t483;
t478 = cos(t482);
t477 = cos(t481);
t476 = cos(t480);
t475 = sin(t482);
t474 = sin(t481);
t473 = sin(t480);
t1 = [t476 * t517 + t477 * t516 + t478 * t515 - m(4) * g(1) + ((-pkin(1) * cos(t495) - pkin(2) * cos(t485) - pkin(3) * t478) * t518 + (-pkin(1) * cos(t494) - pkin(2) * cos(t484) - pkin(3) * t477) * t520 + (-pkin(1) * cos(t493) - pkin(2) * cos(t483) - pkin(3) * t476) * t519) * t521; t473 * t517 + t474 * t516 + t475 * t515 - m(4) * g(2) + ((-pkin(1) * sin(t495) - pkin(2) * sin(t485) - pkin(3) * t475) * t518 + (-pkin(1) * sin(t494) - pkin(2) * sin(t484) - pkin(3) * t474) * t520 + (-pkin(1) * sin(t493) - pkin(2) * sin(t483) - pkin(3) * t473) * t519) * t521; (-0.3e1 * m(2) - 0.3e1 * m(3) - m(4)) * g(3);];
taugX  = t1;
