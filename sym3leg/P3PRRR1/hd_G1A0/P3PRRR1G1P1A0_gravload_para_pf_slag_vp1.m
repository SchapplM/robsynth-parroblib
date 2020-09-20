% Calculate Gravitation load for parallel robot
% P3PRRR1G1P1A0
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G1P1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:48
% EndTime: 2020-03-09 21:14:48
% DurationCPUTime: 0.24s
% Computational Cost: add. (399->75), mult. (425->132), div. (30->5), fcn. (324->18), ass. (0->54)
t526 = pkin(7) + qJ(2,3);
t517 = qJ(3,3) + t526;
t505 = sin(t517);
t508 = cos(t517);
t511 = sin(t526);
t514 = cos(t526);
t549 = 0.1e1 / (t505 * t514 - t511 * t508);
t527 = pkin(7) + qJ(2,2);
t518 = qJ(3,2) + t527;
t506 = sin(t518);
t509 = cos(t518);
t512 = sin(t527);
t515 = cos(t527);
t548 = 0.1e1 / (t506 * t515 - t512 * t509);
t528 = pkin(7) + qJ(2,1);
t519 = qJ(3,1) + t528;
t507 = sin(t519);
t510 = cos(t519);
t513 = sin(t528);
t516 = cos(t528);
t547 = 0.1e1 / (t507 * t516 - t513 * t510);
t533 = 0.1e1 / pkin(2);
t546 = m(3) / pkin(3) * t533;
t529 = legFrame(3,3);
t520 = sin(t529);
t523 = cos(t529);
t499 = -t520 * g(1) + t523 * g(2);
t502 = t523 * g(1) + t520 * g(2);
t481 = -t508 * (rSges(3,1) * t499 - rSges(3,2) * t502) + t505 * (rSges(3,1) * t502 + rSges(3,2) * t499);
t545 = (((-rSges(2,1) * t499 + rSges(2,2) * t502) * t514 + t511 * (rSges(2,1) * t502 + rSges(2,2) * t499)) * m(2) + ((-t499 * t514 + t511 * t502) * pkin(2) + t481) * m(3)) * t549 * t533;
t530 = legFrame(2,3);
t521 = sin(t530);
t524 = cos(t530);
t500 = -t521 * g(1) + t524 * g(2);
t503 = t524 * g(1) + t521 * g(2);
t482 = -t509 * (rSges(3,1) * t500 - rSges(3,2) * t503) + t506 * (rSges(3,1) * t503 + rSges(3,2) * t500);
t544 = (((-rSges(2,1) * t500 + rSges(2,2) * t503) * t515 + t512 * (rSges(2,1) * t503 + rSges(2,2) * t500)) * m(2) + ((-t500 * t515 + t512 * t503) * pkin(2) + t482) * m(3)) * t548 * t533;
t531 = legFrame(1,3);
t522 = sin(t531);
t525 = cos(t531);
t501 = -t522 * g(1) + t525 * g(2);
t504 = t525 * g(1) + t522 * g(2);
t483 = -t510 * (rSges(3,1) * t501 - rSges(3,2) * t504) + t507 * (rSges(3,1) * t504 + rSges(3,2) * t501);
t543 = (((-rSges(2,1) * t501 + rSges(2,2) * t504) * t516 + t513 * (rSges(2,1) * t504 + rSges(2,2) * t501)) * m(2) + ((-t501 * t516 + t513 * t504) * pkin(2) + t483) * m(3)) * t547 * t533;
t542 = t481 * t549 * t546;
t541 = t482 * t548 * t546;
t540 = t483 * t547 * t546;
t539 = t523 * t505 + t520 * t508;
t538 = -t505 * t520 + t523 * t508;
t537 = t524 * t506 + t521 * t509;
t536 = -t506 * t521 + t524 * t509;
t535 = t525 * t507 + t522 * t510;
t534 = -t507 * t522 + t525 * t510;
t1 = [t534 * t543 - (-pkin(2) * (t513 * t522 - t525 * t516) + t534 * pkin(3)) * t540 + t536 * t544 - (-pkin(2) * (t512 * t521 - t524 * t515) + t536 * pkin(3)) * t541 + t538 * t545 - (-pkin(2) * (t511 * t520 - t523 * t514) + t538 * pkin(3)) * t542 - m(4) * g(1); t535 * t543 - (pkin(2) * (t513 * t525 + t522 * t516) + t535 * pkin(3)) * t540 + t537 * t544 - (pkin(2) * (t512 * t524 + t521 * t515) + t537 * pkin(3)) * t541 + t539 * t545 - (pkin(2) * (t511 * t523 + t520 * t514) + t539 * pkin(3)) * t542 - m(4) * g(2); -0.3e1 * g(3) * (m(1) + m(2) + m(3)) - m(4) * g(3);];
taugX  = t1;
