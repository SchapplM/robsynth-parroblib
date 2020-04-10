% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR1G3P1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G3P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G3P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:35
% EndTime: 2020-03-09 21:02:36
% DurationCPUTime: 0.57s
% Computational Cost: add. (198->84), mult. (387->174), div. (126->10), fcn. (447->18), ass. (0->71)
t1503 = legFrame(3,2);
t1488 = sin(t1503);
t1491 = cos(t1503);
t1482 = t1488 * g(1) + t1491 * g(2);
t1485 = t1491 * g(1) - t1488 * g(2);
t1507 = sin(qJ(2,3));
t1513 = cos(qJ(2,3));
t1472 = t1482 * t1513 - t1485 * t1507;
t1494 = 0.1e1 / t1507;
t1545 = t1472 * t1494;
t1504 = legFrame(2,2);
t1489 = sin(t1504);
t1492 = cos(t1504);
t1483 = t1489 * g(1) + t1492 * g(2);
t1486 = t1492 * g(1) - t1489 * g(2);
t1509 = sin(qJ(2,2));
t1515 = cos(qJ(2,2));
t1475 = t1483 * t1515 - t1486 * t1509;
t1495 = 0.1e1 / t1509;
t1544 = t1475 * t1495;
t1505 = legFrame(1,2);
t1490 = sin(t1505);
t1493 = cos(t1505);
t1484 = t1490 * g(1) + t1493 * g(2);
t1487 = t1493 * g(1) - t1490 * g(2);
t1511 = sin(qJ(2,1));
t1517 = cos(qJ(2,1));
t1478 = t1484 * t1517 - t1487 * t1511;
t1496 = 0.1e1 / t1511;
t1543 = t1478 * t1496;
t1542 = t1482 * t1494;
t1541 = t1483 * t1495;
t1540 = t1484 * t1496;
t1512 = cos(qJ(3,3));
t1497 = 0.1e1 / t1512;
t1539 = t1494 * t1497;
t1514 = cos(qJ(3,2));
t1499 = 0.1e1 / t1514;
t1538 = t1495 * t1499;
t1516 = cos(qJ(3,1));
t1501 = 0.1e1 / t1516;
t1537 = t1496 * t1501;
t1536 = t1488 * t1539;
t1535 = t1489 * t1538;
t1534 = t1490 * t1537;
t1533 = t1491 * t1539;
t1532 = t1492 * t1538;
t1531 = t1493 * t1537;
t1506 = sin(qJ(3,3));
t1530 = t1506 * t1539;
t1529 = t1494 / t1512 ^ 2 * t1513;
t1508 = sin(qJ(3,2));
t1528 = t1508 * t1538;
t1527 = t1495 / t1514 ^ 2 * t1515;
t1510 = sin(qJ(3,1));
t1526 = t1510 * t1537;
t1525 = t1496 / t1516 ^ 2 * t1517;
t1524 = t1472 * t1530;
t1523 = t1475 * t1528;
t1522 = t1478 * t1526;
t1521 = t1506 * t1529;
t1520 = t1508 * t1527;
t1519 = t1510 * t1525;
t1518 = 0.1e1 / pkin(2);
t1481 = (g(1) * t1517 + t1511 * g(2)) * t1493 + t1490 * (g(1) * t1511 - g(2) * t1517);
t1480 = (g(1) * t1515 + g(2) * t1509) * t1492 + t1489 * (g(1) * t1509 - g(2) * t1515);
t1479 = (g(1) * t1513 + g(2) * t1507) * t1491 + t1488 * (g(1) * t1507 - g(2) * t1513);
t1477 = t1484 * t1511 + t1487 * t1517;
t1474 = t1483 * t1509 + t1486 * t1515;
t1471 = t1482 * t1507 + t1485 * t1513;
t1 = [(-(t1511 * t1490 + t1493 * t1517) * t1540 - (t1509 * t1489 + t1492 * t1515) * t1541 - (t1507 * t1488 + t1491 * t1513) * t1542) * MDP(1) - g(1) * MDP(12) + ((t1472 * t1533 + t1475 * t1532 + t1478 * t1531) * MDP(3) + (-t1471 * t1533 - t1474 * t1532 - t1477 * t1531) * MDP(4) + (t1491 * t1545 + t1492 * t1544 + t1493 * t1543) * MDP(10) + (-t1491 * t1524 - t1492 * t1523 - t1493 * t1522) * MDP(11)) * t1518; (-(-t1490 * t1517 + t1493 * t1511) * t1540 - (-t1489 * t1515 + t1492 * t1509) * t1541 - (-t1488 * t1513 + t1491 * t1507) * t1542) * MDP(1) - g(2) * MDP(12) + ((-t1472 * t1536 - t1475 * t1535 - t1478 * t1534) * MDP(3) + (t1471 * t1536 + t1474 * t1535 + t1477 * t1534) * MDP(4) + (-t1488 * t1545 - t1489 * t1544 - t1490 * t1543) * MDP(10) + (t1488 * t1524 + t1489 * t1523 + t1490 * t1522) * MDP(11)) * t1518; (-t1482 * t1530 - t1483 * t1528 - t1484 * t1526) * MDP(1) - g(3) * MDP(12) + ((t1472 * t1521 + t1475 * t1520 + t1478 * t1519) * MDP(3) + (-t1471 * t1521 - t1474 * t1520 - t1477 * t1519) * MDP(4) + (t1517 * t1522 + t1501 * (-g(3) * t1516 + t1481 * t1510) + t1515 * t1523 + t1499 * (-g(3) * t1514 + t1480 * t1508) + t1513 * t1524 + t1497 * (-g(3) * t1512 + t1479 * t1506)) * MDP(10) + (-t1510 ^ 2 * t1478 * t1525 + t1501 * (g(3) * t1510 + t1481 * t1516) - t1508 ^ 2 * t1475 * t1527 + t1499 * (g(3) * t1508 + t1480 * t1514) - t1506 ^ 2 * t1472 * t1529 + t1497 * (g(3) * t1506 + t1479 * t1512)) * MDP(11)) * t1518;];
taugX  = t1;
