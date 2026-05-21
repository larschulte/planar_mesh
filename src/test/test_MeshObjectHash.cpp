#include "MeshObject/MeshObject.hpp"

#include <iostream>
#include <unordered_set>

class Test: public MeshObject
{
public:
    const int& get_id() const override
    {
        return id_;
    }

    void set_id(int id)
    {
        id_ = id;
    }

private:
    int id_;
};


int main()
{
    std::shared_ptr<Test> t0 = std::make_shared<Test>();
    std::shared_ptr<Test> t1 = std::make_shared<Test>();
    std::shared_ptr<Test> t2 = std::make_shared<Test>();

    t0->set_id(0);
    t1->set_id(1);
    t2->set_id(2);

    std::unordered_set<std::shared_ptr<Test>, MeshObjectHash> test_set;
    test_set.insert(t0);
    test_set.insert(t1);
    test_set.insert(t2);

    std::cout << "Before delete" << std::endl;
    std::cout << t0->get_id() << std::endl;
    std::cout << t1->get_id() << std::endl;
    std::cout << t2->get_id() << std::endl;

    test_set.erase(t0);
    std::cout << t0->get_id() << std::endl;
    
    test_set.erase(t1);
    std::cout << t1->get_id() << std::endl;

    test_set.erase(t2);
    std::cout << t2->get_id() << std::endl;
    
    return 0;
}